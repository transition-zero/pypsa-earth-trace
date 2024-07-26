import argparse
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
    "'snapshots={{start: {year}-01-01, end: {year}-12-31}}' "
    "'load_options={{scale: {scale}}}' "
    "'run={{name: {iso}/{year}}}' "
    "'cutouts={{cutout: {iso}}}' "
)
NCORES = 1
SNAKEMAKE = (
    f"snakemake --cores {NCORES} "
    f"{TARGET} "
    f"--configfile /mnt/disks/gcs/{BUCKET}/country_configs/config.{{iso}}.yaml "
    f"--config {CONFIG} "
    f"--dry-run"
)

COMMAND = (
    f"python submit_job "
    f"--command {SNAKEMAKE} "
    f"--image {IMAGE} --image-tag {IMAGE_TAG} "
    f"--gcs-bucket-path {BUCKET} "
    f"--configfile {CONFIGFILE} "
    f"--project-id {PROJECT_ID} "
    f"--region {REGION} "
    f"--machine-type {{machine_type}} "
    f"--disk-size-gb {{disk_size_gb}} "
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--weather-year", type=int, required=True)
    args = parser.parse_args()

    configs = os.listdir("../configs")
    iso_codes = [re.search("[A-Z]{2}", config).group() for config in configs]

    year = args.weather_year
    df_demand_scale_factors = pd.read_csv(
        f"../trace_data/demand_scale_factors_{year}.csv",
        keep_default_na=keep_default_na,
        na_values=na_values,
    )

    for iso in iso_codes:
        print(iso)
        machine_type = "n1-standard-4"  # TODO: determine based on iso_code
        disk_size_gb = 32  # TODO: determine based on iso_code
        scale = df_demand_scale_factors.loc[lambda x: x.iso2 == iso, "scale_factor"].item()

        command = COMMAND.format(
            iso=iso,
            year=year,
            scale=scale,
            machine_type=machine_type,
            disk_size_gb=disk_size_gb,
        )
        p = subprocess.Popen(shlex.split(command))
