import argparse
import ast
import json
import logging
import re
import sys
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"
BUCKET = "feo-pypsa-staging"
IMAGE = "europe-west2-docker.pkg.dev/{tz-feo-staging}/feo-pypsa/pypsa-earth-image:latest"

CONFIGFILE_TEMPLATE = r"config.[A-Z][A-Z].yaml"

SNAKEMAKE = "snakemake {target} --configfile {configs_dir}/config.{iso}.yaml {config} {extra}"


def extract_filename_iso(filename):
    if re.fullmatch(CONFIGFILE_TEMPLATE, filename) is None:
        raise ValueError(
            f"Filename '{filename}' does not match expected template '{CONFIGFILE_TEMPLATE}'."
        )
    return re.search(r"[A-Z]{2}", filename).group()


def check_iso_exists(iso, all_iso):
    if iso not in all_iso:
        logger.info(f"Ignoring included ISO code '{iso}' as it is not found in config files.")
    return iso


def available_iso_codes(configs_dir):
    return [extract_filename_iso(f.name) for f in Path(configs_dir).glob(CONFIGFILE_TEMPLATE)]


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
    if next(Path(arg).glob(CONFIGFILE_TEMPLATE), None) is None:
        raise argparse.ArgumentTypeError("No config files found in directory.")
    return Path(arg)


def parse_snakemake_extra_args(arg):
    arg = arg.strip()
    if "--cores" in arg:
        return arg
    return f"--cores 1 {arg}"


def config_dict_to_snakemake(config_dict):
    return "--config " + " ".join([f"'{k}={json.dumps(v)}'" for k, v in config_dict.items()])


class ParseNestedDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        dest, *keys = option_string.lstrip("--").split("-")
        nested_dict = ParseNestedDict.build_nested_dict(keys, values)

        if hasattr(namespace, dest):
            ParseNestedDict.merge_nested_dicts(getattr(namespace, dest), nested_dict)
        else:
            setattr(namespace, dest, nested_dict)

    @staticmethod
    def build_nested_dict(keys, value):
        if len(keys) == 1:
            return {keys[0]: ast.literal_eval(value)}
        else:
            return {keys[0]: ParseNestedDict.build_nested_dict(keys[1:], value)}

    @staticmethod
    def merge_nested_dicts(d1, d2):
        for k, v in d2.items():
            if k in d1 and isinstance(d1[k], dict) and isinstance(v, dict):
                ParseNestedDict.merge_nested_dicts(d1[k], v)
            else:
                d1[k] = v


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
    # path to directory containing Snakemake configfiles
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

    argv = iter(sys.argv[1:])
    args = []
    for arg in argv:
        # Dynamically add arguments starting with --config-
        if arg.startswith("--config-"):
            arg = arg.split("=")[0]  # Handle case where argument specified as foo=bar
            parser.add_argument(arg, type=str, action=ParseNestedDict)
        # Handle --snakemake-extra-args argument special case with single argument e.g.
        # --snakemake-extra-args '--dry-run' breaks argparse unless you add a space before/after
        if arg == "--snakemake-extra-args":
            args.append(arg)
            snakemake_extra_args = next(argv, None)
            args.append(f" {snakemake_extra_args}")
            continue
        args.append(arg)

    args = parser.parse_args(args)
    return args


def main():
    args = parse_args()
    iso_codes = chosen_iso_codes(args.configs_dir, args.iso_include, args.iso_exclude)
    config = config_dict_to_snakemake(args.config) if "config" in args else ""
    print(args.snakemake_extra_args)
    for iso in iso_codes:
        snakemake = SNAKEMAKE.format(
            target=args.target,
            configs_dir=args.configs_dir,
            iso=iso,
            config=config,
            extra=args.snakemake_extra_args,
        )
        print(snakemake)


if __name__ == "__main__":
    main()
