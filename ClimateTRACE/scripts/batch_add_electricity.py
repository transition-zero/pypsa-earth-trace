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

TARGET = "add_electricity"
CONFIGFILE = "./ClimateTRACE/configs/config.{iso}.yaml"
CONFIG = (
    "'snapshots={{start: {year}-01-01, end: {year_}-01-01}}' "
    "'load_options={{scale: {scale}}}' "
    "'run={{name: {iso}/{year}}}'"
)
SNAKEMAKE = (
    f"snakemake "
    f"--cores 2 "
    f"{TARGET} "
    f"--configfile /mnt/disks/gcs/{BUCKET}/country_configs/config.{{iso}}.yaml "
    # f"--configfile {CONFIGFILE} "
    f"--config {CONFIG} "
    "--dry-run"
)

COMMAND = (
    f'python ./ClimateTRACE/scripts/submit_job '
    f'--image {IMAGE} '
    f'--image-tag {IMAGE_TAG} '
    f'--gcs-bucket-path {BUCKET} '
    f'--project-id {PROJECT_ID} '
    f'--region {REGION} '
    f'--machine-type {{machine_type}} '
    f'--disk-size-gb {{disk_size_gb}} '
    f'--configfile {CONFIGFILE} '
    f'--command "{SNAKEMAKE}"'  
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--weather-year", type=int, required=True)
    parser.add_argument("--iso-codes", nargs="+", required=False)
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
        machine_type = "n1-standard-4"  # TODO: determine based on iso_code
        disk_size_gb = 32  # TODO: determine based on iso_code
        scale = df_demand_scale_factors.loc[lambda x: x.iso2 == iso, "scale_factor"].item()

        command = COMMAND.format(
            iso=iso,
            year=year,
            year_=year + 1,
            scale=1.0 if math.isnan(scale) else scale,
            machine_type=machine_type,
            disk_size_gb=disk_size_gb,
        )
        p = subprocess.run(shlex.split(command))
