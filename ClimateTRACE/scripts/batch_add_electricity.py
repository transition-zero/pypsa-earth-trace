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

TARGET = "solve_all_networks"
CONFIGFILE = "./ClimateTRACE/configs/config.{iso}.yaml"
CONFIG = (
    "'snapshots={{start: {year}-01-01, end: {year_}-01-01}}' "
    "'load_options={{scale: {scale}}}' "
    "'run={{name: {iso}/{year}}}' "
    "'scenario={{ll: ['v1.25'], clusters: [10], opts: [1H]}}'"
)
SNAKEMAKE = (
    f"snakemake "
    f"--cores 1 "
    f"{TARGET} "
    f"--configfile {CONFIGFILE} "
    f"--config {CONFIG} "
    "--dry-run"
)

SYMLINK = (
    f"ln -s /mnt/disks/gcs/{BUCKET}/data . && "
    f"ln -s /mnt/disks/gcs/{BUCKET}/ClimateTRACE ."
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
    "--no-spot "
    f"--configfile {CONFIGFILE} "
    f"--command {COMMAND}"
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--weather-year", type=int, required=True)
    parser.add_argument("--iso-codes", nargs="+", required=False)
    parser.add_argument("--local", action="store_true", default=False)
    parser.add_argument("--dry-run", action="store_true", default=False)
    args = parser.parse_args()

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
    for iso in iso_codes:
        print(f"this is the iso: {iso}")
        machine_type = "n1-standard-8"  # TODO: scale based on iso_code
        disk_size_gb = 64  # TODO: scale based on iso_code
        try:
            scale = df_demand_scale_factors.loc[
                lambda x: x.iso2 == iso, "scale_factor"
            ].item()
        except ValueError:
            scale = 1.0
        submit = SUBMIT.format(
            iso=iso,
            year=year,
            year_=year + 1,
            scale=1.0 if math.isnan(scale) else scale,
            machine_type=machine_type,
            disk_size_gb=disk_size_gb,
        )

        snakemake = re.search(r"snakemake .*", submit).group(0)
        if args.dry_run:
            snakemake += " --dry-run"
            submit = re.sub(r"snakemake .*", snakemake, submit)

        if args.local:
            print(f"Running locally with command: {snakemake}")
            run_args = shlex.split(snakemake)
        else:
            command = re.search(r"--command (.*)", submit).group(1)
            # NOTE: shlex.split doesn't work nicely on nested shell statements without very specific
            # quote escaping, if at all. It's cleaner to extract the nested command before splitting
            # and handle it cleanly separetely in the submit_job script.
            run_args = shlex.split(submit.replace(command, "")) + [command]

        p = subprocess.run(run_args)
