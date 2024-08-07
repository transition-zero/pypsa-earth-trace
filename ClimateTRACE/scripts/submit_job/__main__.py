# -*- coding: utf-8 -*-
# mypy: ignore-errors

import os
import re
import shlex
import yaml
from datetime import datetime

import click
from batch_job import create_container_job, wait_for_jobs_to_succeed
from gcs_file_utils import upload_file_to_bucket


@click.group()
def cli():
    pass


class YamlParam(click.ParamType):
    def convert(self, value, param, ctx):
        try:
            return yaml.load(value, Loader=yaml.BaseLoader)
        except yaml.YAMLError as e:
            self.fail(f"Could not parse YAML: {e}")


@cli.command("submit_job")
@click.option("--project-id")
@click.option("--region")
@click.option("--image")
@click.option("--image-tag")
@click.option("--command")
@click.option("--max-retries", default=0)
@click.option("--max-duration", default="27000s")
@click.option("--task-environments", "-e", multiple=True, type=YamlParam())
@click.option("--parallelism", default=1, type=int)
@click.option("--machine-type", type=str, default="n1-standard-4")
@click.option("--no-spot", is_flag=True)
@click.option(
    "--disk-size-gb",
    type=int,
    default=None,
    help=(
        "Boot disk size in GB. Batch will calculate the boot disk size based on "
        "source image and task requirements if you do not speicify the size."
    ),
)
@click.option("--gcs-bucket-path", type=str)
@click.option("--configfiles", "-f", multiple=True, type=str)
@click.option("--wait-for-completion", is_flag=True)
def submit_job(
    project_id: str,
    region: str,
    image: str,
    image_tag: str,
    command: str,
    max_retries: int,
    max_duration: str,
    task_environments: list[dict[str, str]],
    parallelism: int,
    machine_type: str,
    no_spot: bool,
    disk_size_gb: int | None,
    gcs_bucket_path: str | None,
    configfiles: list[str] | None,
    wait_for_completion: bool,
):
    print("starting submit job")
    for configfile in configfiles:
        upload_file_to_bucket(
            bucket_name=gcs_bucket_path,
            blob_name=f"ClimateTRACE/configs/{os.path.basename(configfile)}",
            local_file_name=configfile,
            content_type="application/x-yaml",
        )

    snakemake = re.search(r"snakemake .*", command).group(0)
    iso = "" if not (match := re.search(r"config.([A-Z]{2}).yaml", snakemake)) else match.group(1)
    job_id = "-".join(
        [
            f"{iso.lower()}",
            snakemake.lower()
            .split(" ")[3]
            .split("/")[-1]
            .replace(".", "-")
            .replace("_", "-")
            .strip("-"),
            f"{datetime.now().strftime('%Y%m%d%H%M%S')}",
        ]
    ).lstrip("-")
    print(f"job_id: {job_id}")

    if re.match("/bin/bash -c ", command):
        commands = re.split(r"^(/bin/bash) (-c) ", command)[1:]
    else:
        commands = shlex.split(command)
    print(f"commands: {commands}")

    job = create_container_job(
        project_id=project_id,
        region=region,
        image=image,
        image_tag=image_tag,
        commands=commands,
        max_retries=max_retries,
        max_duration=max_duration,
        task_environments=task_environments,
        parallelism=parallelism,
        machine_type=machine_type,
        spot=(not no_spot),
        disk_size_gb=disk_size_gb,
        gcs_bucket_path=gcs_bucket_path,
        job_id=job_id,
    )
    print(f"Starting job '{job.name}'")
    if wait_for_completion:
        wait_for_jobs_to_succeed([job])


if __name__ == "__main__":
    submit_job()
