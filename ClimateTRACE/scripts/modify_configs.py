# -*- coding: utf-8 -*-
import argparse
from pathlib import Path

from ruamel.yaml import YAML


def load_yaml(file_path: str):
    yaml = YAML()
    yaml.preserve_quotes = True
    with open(file_path, "r") as file:
        return yaml, yaml.load(file)


def save_yaml(file_path: str, yaml, data: YAML):
    with open(file_path, "w") as file:
        yaml.dump(data, file)


def change_snapshots(file_path: str, year_start: int):
    yaml, data = load_yaml(file_path)
    year_start_str = f"{year_start}-01-01"
    year_end_str = f"{year_start + 1}-01-01"
    if "snapshots" in data and isinstance(data["snapshots"], dict):
        snapshot = data["snapshots"]
        snapshot["start"] = year_start_str
        snapshot["end"] = year_end_str
    save_yaml(file_path, yaml, data)


def change_load_scale(file_path, new_load_scale):
    yaml, data = load_yaml(file_path)
    if "load_options" in data and "scale" and isinstance(data["load_options"], dict):
        data["load_options"]["scale"] = new_load_scale
    save_yaml(file_path, yaml, data)


def change_enable_rules(file_path):
    yaml, data = load_yaml(file_path)
    if "enable" in data and isinstance(data["enable"], dict):
        data["enable"]["retrieve_databundle"] = False
        data["enable"]["retrieve_cost_data"] = False
        data["enable"]["download_osm_data"] = False
        data["enable"]["build_natura_raster"] = False
        data["enable"]["build_cutout"] = False
    save_yaml(file_path, yaml, data)


def modify_config(
    file_path: str,
    year_start: int = None,
    load_scale: float = None,
    enable_rules: bool = False,
):
    if year_start is not None:
        change_snapshots(
            file_path=file_path,
            year_start=year_start,
        )
    if load_scale is not None:
        change_load_scale(file_path=file_path, new_load_scale=load_scale)
    if enable_rules:
        change_enable_rules(file_path=file_path)


def main():
    parser = argparse.ArgumentParser(description="Update YAML configuration file.")
    parser.add_argument(
        "config_path",
        help="Directory containing YAML configuration files, or path to a single config file.",
        type=str,
    )
    parser.add_argument("--year_start", help="New start date for snapshots.", type=int)
    parser.add_argument(
        "--load_scale",
        type=float,
        help="Scales load up or down",
    )
    parser.add_argument(
        "--enable_rules", action="store_true", help="Change enable rules."
    )

    args = parser.parse_args()

    config_path = Path(args.config_path)
    if config_path.is_file():
        print(f"Modifying {config_path}")
        modify_config(
            file_path=config_path,
            year_start=args.year_start,
            load_scale=args.load_scale,
            enable_rules=args.enable_rules,
        )
    elif config_path.is_dir():
        for file_path in config_path.glob("config.*.yaml"):
            print(f"Modifying {file_path}")
            modify_config(
                file_path=file_path,
                year_start=args.year_start,
                load_scale=args.load_scale,
                enable_rules=args.enable_rules,
            )
    else:
        raise FileNotFoundError(f"File or directory not found: {config_path}")


if __name__ == "__main__":
    main()
