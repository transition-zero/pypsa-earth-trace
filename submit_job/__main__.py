# -*- coding: utf-8 -*-
# mypy: ignore-errors

import os

import click
from batch_job import create_container_job, wait_for_jobs_to_succeed
from gcs_file_utils import upload_file_to_bucket


@click.group()
def cli():
    pass


@cli.command("submit_job")
@click.option("--project-id", default="tz-feo-staging")
@click.option("--region", default="europe-west2")
@click.option("--image")
@click.option("--image-tag")
@click.option("--command")
@click.option("--max-retries", default=0)
@click.option("--max-duration", default="3600s")
@click.option("--task-count", default=1, type=int)
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
@click.option(
    "--cpu-milli",
    type=int,
    default=None,
    help="CPU request per task. Defaults to all CPUs available on selected VM.",
)
@click.option(
    "--memory-mb",
    type=int,
    default=None,
    help="Memory request per task. Defaults to all RAM available on selected VM.",
)
@click.option("--gcs-bucket-path", type=str, default="feo-pypsa-staging")
@click.option("--config-file", type=str, default="config.default.yaml")
@click.option("--wait-for-completion", is_flag=True)
def submit_job(
    project_id: str,
    region: str,
    image: str,
    image_tag: str,
    command: str,
    max_retries: int,
    max_duration: str,
    task_count: int,
    parallelism: int,
    machine_type: str,
    no_spot: bool,
    disk_size_gb: int | None,
    cpu_milli: int | None,
    memory_mb: int | None,
    gcs_bucket_path: str | None,
    config_file: str | None,
    wait_for_completion: bool,
):
    print("starting submit job")
    upload_file_to_bucket(
        bucket_name=gcs_bucket_path,
        blob_name=f"country_configs/{os.path.basename(config_file)}",
        local_file_name=config_file,
        content_type="application/x-yaml",
    )
    job = create_container_job(
        project_id=project_id,
        region=region,
        image=image,
        image_tag=image_tag,
        command=command,
        max_retries=max_retries,
        max_duration=max_duration,
        task_count=task_count,
        parallelism=parallelism,
        machine_type=machine_type,
        spot=(not no_spot),
        disk_size_gb=disk_size_gb,
        cpu_milli_per_task=cpu_milli,
        memory_mb_per_task=memory_mb,
        gcs_bucket_path=gcs_bucket_path,
    )
    if wait_for_completion:
        wait_for_jobs_to_succeed([job])
    print("job submitted")


if __name__ == "__main__":
    submit_job()
