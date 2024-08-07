import argparse
import math
import os
import re
import shlex
import subprocess

import pandas as pd
from pandas._libs.parsers import STR_NA_VALUES

try:
    STR_NA_VALUES.remove("NA")
    keep_default_na = False
    na_values = STR_NA_VALUES
except KeyError:
    keep_default_na = True
    na_values = None

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"
BUCKET = "feo-pypsa-staging"
IMAGE = "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image"
IMAGE_TAG = "latest"

CONFIGFILE = "./ClimateTRACE/configs/config.$ISO_COUNTRY_CODE.yaml"
CONFIG = (
    "'snapshots={{start: {year}-01-01, end: {year_}-01-01}}' "
    "'load_options={{scale: $ISO_DEMAND_SCALE}}' "
    "'run={{name: $ISO_COUNTRY_CODE/{year}}}' "
    "'scenario={{simpl: [\"\"], ll: ['v1.25'], clusters: ['1'], opts: ['1H']}}' "
)
SNAKEMAKE = (
    f"snakemake "
    f"--cores 1 "
    "{target} "
    f"--configfile {CONFIGFILE} "
    f"--config {CONFIG} "
    # "--rerun-triggers mtime "
)

SYMLINK = (
    f"ln -s /mnt/disks/gcs/{BUCKET}/data . && "
    f"ln -s /mnt/disks/gcs/{BUCKET}/ClimateTRACE . "
)
COMMAND = f"/bin/bash -c {SYMLINK} && {SNAKEMAKE}"
SUBMIT = (
    f"python ./ClimateTRACE/scripts/submit_job "
    f"--image {IMAGE} "
    f"--image-tag {IMAGE_TAG} "
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=str, required=True)
    parser.add_argument("--weather-year", type=int, required=True)
    parser.add_argument("--iso-codes", nargs="+", required=False)
    parser.add_argument("--local", action="store_true", default=False)
    parser.add_argument("--dry-run", action="store_true", default=False)
    parser.add_argument("--machine-type", type=str, default="n1-standard-8")
    parser.add_argument("--disk-size-gb", type=int, default=64)
    parser.add_argument("--parallelism", type=int, default=1)
    # TODO: enable running each country as a separate batch job (i.e. like before)
    # parser.add_argument("--batch-mode", type=str, choices=["tasks", "jobs"], default="jobs")
    args = parser.parse_args()

    if args.dry_run:
        SNAKEMAKE += "--dry-run "
        SUBMIT += "--dry-run "

    configs = os.listdir("./ClimateTRACE/configs")
    iso_codes = (
        args.iso_codes
        if args.iso_codes
        else [re.search("[A-Z]{2}", config).group() for config in configs]
    )

    year = args.weather_year

    df_demand_scale_factors = pd.read_csv(
        f"./ClimateTRACE/trace_data/demand_scale_factors_{year}.csv",
        keep_default_na=keep_default_na,
        na_values=na_values,
    )
    print(iso_codes)

    task_environments = []
    for iso in iso_codes:
        print(f"this is the iso: {iso}")
        try:
            scale = df_demand_scale_factors.loc[lambda x: x.iso2 == iso, "scale_factor"].item()
        except ValueError:
            scale = 1.0
        scale = 1.0 if math.isnan(scale) else scale

        if args.local:
            snakemake = (
                SNAKEMAKE.format(target=args.target, year=year, year_=year + 1)
                .replace("$ISO_COUNTRY_CODE", iso)
                .replace("$ISO_DEMAND_SCALE", str(scale))
            )
            print(f"Running locally with command: {snakemake}")
            run_args = shlex.split(snakemake)
            p = subprocess.run(run_args)
        else:
            task_environments.append(f"-e '{{ISO_COUNTRY_CODE: {iso}, ISO_DEMAND_SCALE: {scale}}}'")

    if not args.local:
        submit = SUBMIT.format(
            target=args.target,
            task_environments=" ".join(task_environments),
            configfiles=" ".join(
                [
                    f"-f {CONFIGFILE.replace('$ISO_COUNTRY_CODE', iso)}"
                    for iso in iso_codes
                ]
            ),
            year=year,
            year_=year + 1,
            machine_type=args.machine_type,
            disk_size_gb=args.disk_size_gb,
            parallelism=args.parallelism,
        )

        command = re.search(r"--command (.*)", submit).group(1)
        # NOTE: shlex.split doesn't work nicely on nested shell statements without very specific
        # quote escaping, if at all. It's cleaner to extract the nested command before splitting
        # and handle it cleanly separetely in the submit_job script.
        run_args = shlex.split(submit.replace(command, "")) + [command]

        p = subprocess.run(run_args)
