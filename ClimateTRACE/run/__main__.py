import argparse
import math
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime

import gcp_utils
import pandas as pd

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"
BUCKET = "feo-pypsa-staging"
IMAGE_URI = "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image:latest"

CONFIGFILE = "ClimateTRACE/configs/config.{iso}.yaml"
# NOTE: the order and the type of quotes is important downstream for the shell to be able to
# correctly perform variable substitution within the json strings.
CONFIG = (
    '"snapshots={{start: {year}-01-01, end: {year_}-01-01}}" '
    '"load_options={{scale: {scale}}}" '
    '"run={{name: {iso}/{year}}}" '
    "\"scenario={{simpl: [''], ll: [v1.25], clusters: [1], opts: [1H]}}\" "
)
SNAKEMAKE = (
    f"snakemake "
    "{target} "
    f"--configfile {CONFIGFILE} "
    f"--config {CONFIG} "
    "{snakemake_extra_args} "
)
# TODO: get rid of the symlink and mount gcs volumes directly in container workdir
SYMLINK = f"ln -s /mnt/disks/gcs/{BUCKET}/data . && ln -s /mnt/disks/gcs/{BUCKET}/ClimateTRACE ."


def submit_job(
    job_id: str,
    command: str,
    task_environments: list[dict[str, str]],
    parallelism: int,
    machine_type: str,
    no_spot: bool,
    disk_size_gb: int,
    configfiles: list[str],
):
    if re.match("/bin/bash -c ", command):
        commands = re.split(r"^(/bin/bash) (-c) ", command)[1:]
    else:
        commands = shlex.split(command)

    print(f"job_id: {job_id}")
    print(f"commands: {commands}")
    print(f"task_environments: {task_environments}")
    print(f"parallelism: {parallelism}")

    for configfile in configfiles:
        gcp_utils.upload_file_to_bucket(
            bucket_name=BUCKET,
            blob_name=configfile,
            local_file_name=configfile,
            content_type="application/x-yaml",
        )
    job = gcp_utils.create_container_job(
        project_id=PROJECT_ID,
        region=REGION,
        image_uri=IMAGE_URI,
        commands=commands,
        task_environments=task_environments,
        parallelism=parallelism,
        machine_type=machine_type,
        spot=(not no_spot),
        disk_size_gb=disk_size_gb,
        gcs_bucket_path=BUCKET,
        job_id=job_id,
    )
    print(f"Starting job {job.name}{os.linesep}")
    return job


def snakemake_job_id(snakemake):
    iso = "" if not (match := re.search(r"config.([A-Z]{2}).yaml", snakemake)) else match.group(1)
    job_id = "-".join(
        [
            f"{iso.lower()}",
            snakemake.lower()
            .split(" ")[1]
            .split("/")[-1]
            .replace(".", "-")
            .replace("_", "-")
            .strip("-"),
            f"{datetime.now().strftime('%Y%m%d%H%M%S')}",
        ]
    ).lstrip("-")
    return job_id


def demand_scale_factor(iso, df_demand_scale_factors):
    try:
        scale = df_demand_scale_factors.loc[lambda x: x.iso2 == iso, "scale_factor"].item()
    except ValueError:
        scale = 1.0
    return 1.0 if math.isnan(scale) else scale


def chosen_iso_codes(iso_include=None, iso_exclude=None):
    all_iso = [
        re.search("[A-Z]{2}", config).group() for config in os.listdir("./ClimateTRACE/configs")
    ]
    if iso_include is not None:
        if "all" in iso_include:
            return all_iso
        return [iso for iso in iso_include if iso in all_iso]
    if iso_exclude is not None:
        return [iso for iso in all_iso if iso not in iso_exclude]


def parse_args():
    parser = argparse.ArgumentParser()

    smk = parser.add_argument_group("snakemake arguments")
    iso = smk.add_mutually_exclusive_group(required=True)
    smk.add_argument(
        "--target", type=str, required=True, help="Snakemake target rule/file name to run."
    )
    smk.add_argument(
        "--snakemake-extra-args",
        type=lambda arg_: f"--cores 1 {arg}" if "--cores" not in (arg := arg_.strip()) else arg,
        default="--cores 1",
        required=False,
        help="Extra arguments to pass to Snakemake.",
    )
    iso.add_argument(
        "--iso-include",
        type=str,
        nargs="+",
        help="Country ISO codes to INCLUDE from config files, or 'all'.",
    )
    iso.add_argument(
        "--iso-exclude",
        type=str,
        nargs="+",
        help="Country ISO codes to EXCLUDE from config files",
    )
    smk.add_argument("--weather-year", type=int, required=True, help="Model weather year to use.")

    job = parser.add_argument_group("job configuration")
    local_or_batch = job.add_mutually_exclusive_group(required=True)
    local_or_batch.add_argument(
        "--local", action="store_true", default=False, help="Run snakemake workflow locally."
    )
    local_or_batch.add_argument(
        "--batch-mode",
        type=str,
        choices=["tasks", "jobs"],
        help="Run each task as a separate job or all tasks as a single job.",
    )

    gcp = parser.add_argument_group("gcloud arguments")
    gcp.add_argument("--machine-type", type=str, default="n1-standard-8", help="VM machine type.")
    gcp.add_argument("--disk-size-gb", type=int, default=64, help="VM boot disk size in GB.")
    gcp.add_argument(
        "--task-parallelism",
        type=int,
        default=1,
        help="Number of parallel tasks to run in batch mode 'tasks' (ignored in 'jobs' mode).",
    )

    # Handle --snakemake-extra-args argument special case with single argument e.g.
    # --snakemake-extra-args '--dry-run' breaks, but --snakemake-extra-args '--dry-run ' works
    argv = iter(sys.argv[1:])
    args = []
    for arg in argv:
        if arg == "--snakemake-extra-args":
            args.append(arg)
            snakemake_extra_args = next(argv, None)
            args.append(f"{snakemake_extra_args} ")
            continue
        args.append(arg)

    return parser.parse_args(args)


if __name__ == "__main__":
    if os.path.split(os.getcwd())[-1] != "pypsa-earth-trace":
        raise RuntimeError("Please run this script from the pypsa-earth-trace directory")

    args = parse_args()

    df_demand_scale_factors = pd.read_csv(
        f"ClimateTRACE/trace_data/demand_scale_factors_{args.weather_year}.csv",
        keep_default_na=False,
        na_values=list(filter(lambda x: x != "NA", pd._libs.parsers.STR_NA_VALUES)),
    )
    iso_codes = chosen_iso_codes(args.iso_include, args.iso_exclude)
    print(f"Running {args.target} in {args.weather_year} for {iso_codes}")

    task_environments = []
    for iso in iso_codes:
        scale = demand_scale_factor(iso, df_demand_scale_factors)
        snakemake = SNAKEMAKE.format(
            target=args.target,
            iso=iso,
            scale=scale,
            year=args.weather_year,
            year_=args.weather_year + 1,
            snakemake_extra_args=args.snakemake_extra_args,
        )
        if args.local:
            print("\n", f"Running locally with command: {snakemake}", "\n")
            subprocess.run(shlex.split(snakemake))
        elif args.batch_mode == "jobs":
            submit_job(
                job_id=snakemake_job_id(snakemake),
                command=f"/bin/bash -c {SYMLINK} && {snakemake}",
                task_environments=None,
                parallelism=1,
                machine_type=args.machine_type,
                no_spot=True,
                disk_size_gb=args.disk_size_gb,
                configfiles=[CONFIGFILE.format(iso=iso)],
            )
        elif args.batch_mode == "tasks":
            task_environments.append({"ISO_COUNTRY_CODE": iso, "ISO_DEMAND_SCALE": str(scale)})

    if not args.local and args.batch_mode == "tasks":
        snakemake = SNAKEMAKE.format(
            target=args.target,
            iso="$ISO_COUNTRY_CODE",
            scale="$ISO_DEMAND_SCALE",
            year=args.weather_year,
            year_=args.weather_year + 1,
            snakemake_extra_args=args.snakemake_extra_args,
        )
        submit_job(
            job_id=snakemake_job_id(snakemake),
            command=f"/bin/bash -c {SYMLINK} && {snakemake}",
            task_environments=task_environments,
            parallelism=args.task_parallelism,
            machine_type=args.machine_type,
            no_spot=True,
            disk_size_gb=args.disk_size_gb,
            configfiles=[CONFIGFILE.format(iso=iso) for iso in iso_codes],
        )
