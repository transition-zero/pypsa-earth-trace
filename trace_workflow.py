import argparse
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"
BUCKET = "feo-pypsa-staging"
IMAGE = "europe-west2-docker.pkg.dev/{tz-feo-staging}/feo-pypsa/pypsa-earth-image:latest"

SNAKEMAKE = "snakemake {target} --configfile {configs_dir}/config.{iso}.yaml --config {{config}}"


def extract_filename_iso(filename):
    return re.search(r"config\.([A-Z]{2})\.yaml", filename).group(1)


def check_iso_exists(iso, all_iso):
    if iso not in all_iso:
        logger.info(f"Ignoring included ISO code '{iso}' as it is not found in config files.")
    return iso


def available_iso_codes(configs_dir):
    return [extract_filename_iso(f.as_posix()) for f in Path(configs_dir).glob("config.*.yaml")]


def chosen_iso_codes(configs_dir, iso_include=None, iso_exclude=None):
    all_iso = available_iso_codes(configs_dir)
    if iso_include is not None:
        if "all" in iso_include:
            return all_iso
        return [iso for iso in iso_include if check_iso_exists(iso, all_iso) in all_iso]
    if iso_exclude is not None:
        return [iso for iso in all_iso if iso not in iso_exclude]


def validate_iso2(arg):
    if re.fullmatch(r"[A-Z]{2}", arg) is None:
        raise argparse.ArgumentTypeError("Argument must be a 2 letter ISO code.")
    return arg


def validate_iso2_all(arg):
    if arg == "all":
        return arg
    try:
        arg = validate_iso2(arg)
        return arg
    except argparse.ArgumentTypeError:
        raise argparse.ArgumentTypeError("Argument must be a 2 letter ISO code or 'all'.")


def validate_config_dir(arg):
    if not Path(arg).is_dir():
        raise argparse.ArgumentTypeError("Argument must be a valid directory path.")
    if next(Path(arg).glob("config.*.yaml"), None) is None:
        raise argparse.ArgumentTypeError("No config files found in directory.")
    return Path(arg)


def parse_snakemake_extra_args(arg):
    if "--cores" in arg:
        return arg
    return f"--cores 1 {arg}"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ClimateTRACE snakemake workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    required_args = parser.add_argument_group("REQUIRED ARGUMENTS")
    optional_args = parser.add_argument_group("OPTIONAL ARGUMENTS")
    iso_args = required_args.add_mutually_exclusive_group(required=True)

    # snakemake target rule/file
    required_args.add_argument(
        "--target", type=str, required=True, help="Snakemake target rule/file name."
    )
    # path to directory containing config files
    required_args.add_argument(
        "--configs-dir",
        type=validate_config_dir,
        required=True,
        help="Path to Snakemake configfiles.",
    )
    # country iso codes to include
    iso_args.add_argument(
        "--iso-include",
        type=validate_iso2_all,
        nargs="+",
        help="Country ISO codes to INCLUDE from config files in --configs-dir, or 'all'.",
    )
    # country iso codes to exclude
    iso_args.add_argument(
        "--iso-exclude",
        type=validate_iso2,
        nargs="+",
        help="Country ISO codes to EXCLUDE from config files in --configs-dir",
    )
    # run workflow on Google Cloud Platform
    optional_args.add_argument(
        "--gcp",
        action="store_true",
        default=False,
        help="Submit workflow to Google Cloud Platform Batch service.",
    )
    # optional extra arguments passed to snakemake command
    optional_args.add_argument(
        "--snakemake-extra-args",
        type=parse_snakemake_extra_args,
        default="--cores 1",
        help="Extra arguments to snakemake as a string, e.g. '--cores 2 --rerun-triggers mtime'.",
    )
    args, unknown = parser.parse_known_args()
    if unknown:
        logger.warning(f"Unknown arguments: {unknown}")
    return args


def main():
    args = parse_args()
    iso_codes = chosen_iso_codes(args.configs_dir, args.iso_include, args.iso_exclude)
    snakemake = SNAKEMAKE + (f" {args.snakemake_extra_args}" if args.snakemake_extra_args else "")

    for iso in iso_codes:
        _snakemake = snakemake.format(
            target=args.target, n=1, configs_dir=args.configs_dir, iso=iso
        )
        print(_snakemake)


if __name__ == "__main__":
    main()
