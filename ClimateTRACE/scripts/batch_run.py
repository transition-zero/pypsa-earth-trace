import argparse
import math
import os
import re
import shlex
import subprocess
import sys

import pandas as pd

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"
BUCKET = "feo-pypsa-staging"
IMAGE_URI = "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image:latest"

CONFIGFILE = "./ClimateTRACE/configs/config.{iso}.yaml"
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

SYMLINK = f"ln -s /mnt/disks/gcs/{BUCKET}/data . && ln -s /mnt/disks/gcs/{BUCKET}/ClimateTRACE ."
COMMAND = f"/bin/bash -c {SYMLINK} && {SNAKEMAKE}"
SUBMIT = (
    f"python ./ClimateTRACE/scripts/submit_job "
    f"--image-uri {IMAGE_URI} "
    f"--gcs-bucket-path {BUCKET} "
    f"--project-id {PROJECT_ID} "
    f"--region {REGION} "
    "--machine-type {machine_type} "
    "--disk-size-gb {disk_size_gb} "
    "--parallelism {parallelism} "
    "--no-spot "
    "{task_environments} "
    "{configfiles} "
    f"--command {COMMAND} "
)


def submit_job(
    *,
    target,
    snakemake_extra_args,
    iso,
    scale,
    year,
    configfiles=None,
    task_environments=None,
    machine_type=None,
    disk_size_gb=None,
    task_parallelism=None,
    local=False,
):
    submit = SUBMIT.format(
        target=target,
        iso=iso,
        scale=scale,
        year=year,
        year_=year + 1,
        snakemake_extra_args=snakemake_extra_args,
        task_environments=task_environments,
        configfiles=configfiles,
        machine_type=machine_type,
        disk_size_gb=disk_size_gb,
        parallelism=task_parallelism,
    )
    if local:
        snakemake = re.search(r"snakemake .*$", submit).group(0)
        print("\n", f"Running locally with command: {snakemake}", "\n")
        run_args = shlex.split(snakemake)
    else:
        command = re.search(r"--command (.*)$", submit).group(1)
        # NOTE: shlex.split doesn't work nicely on nested shell statements without very specific
        # quote escaping, if at all. It's cleaner to extract the nested command before splitting
        # and handle it separetely in the submit_job script.
        run_args = shlex.split(submit.replace(command, "")) + [command]
    p = subprocess.run(run_args)
    return p.returncode


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
    args = parse_args()

    df_demand_scale_factors = pd.read_csv(
        f"./ClimateTRACE/trace_data/demand_scale_factors_{args.weather_year}.csv",
        keep_default_na=False,
        na_values=list(filter(lambda x: x != "NA", pd._libs.parsers.STR_NA_VALUES)),
    )

    iso_codes = chosen_iso_codes(args.iso_include, args.iso_exclude)
    print(f"Running {args.target} in {args.weather_year} for {iso_codes}")

    task_environments = ""
    for iso in iso_codes:
        scale = demand_scale_factor(iso, df_demand_scale_factors)

        if args.local:
            submit_job(
                target=args.target,
                snakemake_extra_args=args.snakemake_extra_args,
                iso=iso,
                scale=scale,
                year=args.weather_year,
                local=True,
            )
        elif args.batch_mode == "jobs":
            submit_job(
                target=args.target,
                snakemake_extra_args=args.snakemake_extra_args,
                iso=iso,
                scale=scale,
                year=args.weather_year,
                task_environments=task_environments,
                configfiles=f"-f {CONFIGFILE.format(iso=iso)}",
                machine_type=args.machine_type,
                disk_size_gb=args.disk_size_gb,
                task_parallelism=1,
                local=False,
            )
        else:
            scale = demand_scale_factor(iso, df_demand_scale_factors)
            task_environments += f"-e '{{ISO_COUNTRY_CODE: {iso}, ISO_DEMAND_SCALE: {scale}}}' "

    if not args.local and args.batch_mode == "tasks":
        submit_job(
            target=args.target,
            snakemake_extra_args=args.snakemake_extra_args,
            iso="$ISO_COUNTRY_CODE",
            scale="$ISO_DEMAND_SCALE",
            year=args.weather_year,
            task_environments=task_environments,
            configfiles=" ".join([f"-f {CONFIGFILE.format(iso=iso)}" for iso in iso_codes]),
            machine_type=args.machine_type,
            disk_size_gb=args.disk_size_gb,
            task_parallelism=args.task_parallelism,
            local=False,
        )
