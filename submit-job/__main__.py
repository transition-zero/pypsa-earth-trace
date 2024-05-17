# mypy: ignore-errors


import click
import batch_job


@click.group()
def cli():
    pass


@cli.command("submit-job")
@click.option("--project-id", default="tz-feo-staging")
@click.option("--region", default="europe-west2")
@click.option("--image", default="gcr.io/tz-ml-dev/eo-models-cpu")
@click.option("--image-tag")
@click.option("--command")
@click.option(
    "--machine-type",
    type=str,
    default="n1-standard-4",
)
@click.option("--max-retries", default=0)
@click.option("--max-duration", default="10000s")
@click.option("--gpu-type", type=str, default=None)
@click.option("--disk-size-gb", default=100, type=int)
@click.option("--task-count", default=1, type=int)
@click.option("--parallelism", default=1, type=int)
@click.option("--no-spot", is_flag=True)
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
@click.option("--wait-for-completion", is_flag=True)
@click.option("--gcs-bucket-path", type=str, default="feo-pypsa-staging")
@click.option("--config-file", type=str, default="config.default.yaml")
def submit_job(
    project_id: str,
    region: str,
    image: str,
    image_tag: str,
    command: str,
    max_retries: int,
    max_duration: str,
    machine_type: str,
    gpu_type: str | None,
    disk_size_gb: int,
    task_count: int,
    parallelism: int,
    no_spot: bool,
    cpu_milli: int | None,
    memory_mb: int | None,
    wait_for_completion: bool,
    gcs_bucket_path: str | None,
    config_file: str | None,
):
    print("starting submit job")
    job = batch_job.create_container_job(
        project_id=project_id,
        region=region,
        image=image,
        image_tag=image_tag,
        command=command,
        max_retries=max_retries,
        max_duration=max_duration,
        machine_type=machine_type,
        gpu_type=gpu_type,
        disk_size_gb=disk_size_gb,
        parallelism=parallelism,
        task_count=task_count,
        spot=(not no_spot),
        cpu_milli_per_task=cpu_milli,
        memory_mb_per_task=memory_mb,
        gcs_bucket_path=gcs_bucket_path,
        config_file=config_file,
    )
    if wait_for_completion:
        batch_job.wait_for_jobs_to_succeed([job])
    print("job submitted")


if __name__ == "__main__":
    submit_job()
