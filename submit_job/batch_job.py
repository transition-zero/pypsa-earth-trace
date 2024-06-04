import os
from collections.abc import Sequence
from time import sleep
from datetime import datetime

from google.cloud import batch_v1
from loguru import logger

from gcs_file_utils import upload_file_to_bucket

# TODO: can we get these from a GCP API?
MACHINE_TYPE_TO_CPU: dict[str, int] = {
    "n1-standard-1": 1_000,
    "n1-standard-2": 2_000,
    "n1-standard-4": 4_000,
    "n1-standard-8": 8_000,
    "n1-standard-16": 16_000,
    "n1-standard-32": 32_000,
    "n1-highmem-1": 1_000,
    "n1-highmem-2": 2_000,
    "n1-highmem-4": 4_000,
    "n1-highmem-8": 8_000,
    "n1-highmem-16": 16_000,
    "n1-highmem-32": 32_000,
    "a2-highgpu-1g": 12_000,
    "n1-highcpu-64": 64_000,
}

MACHINE_TYPE_TO_RAM: dict[str, int] = {
    "n1-standard-1": 3_750,
    "n1-standard-2": 7_500,
    "n1-standard-4": 15_000,
    "n1-standard-8": 30_000,
    "n1-standard-16": 60_000,
    "n1-highmem-1": 6_500,
    "n1-highmem-2": 13_000,
    "n1-highmem-4": 26_000,
    "n1-highmem-8": 52_000,
    "n1-highmem-16": 104_000,
    "n1-highmem-32": 208_000,
    "a2-highgpu-1g": 85_000,
    "n1-highcpu-64": 57_600,
}

BATCH_STATE = batch_v1.JobStatus.State


def wait_for_jobs_to_succeed(batch_jobs: Sequence[batch_v1.Job]):
    """
    Block until all Batch jobs have succeeded. Raise exception if any job fails.

    Args:
        batch_jobs (Sequence[batch_v1.Job]): The Batch jobs to wait for.
    """
    all_jobs_succeeded = False
    while all_jobs_succeeded is False:
        with batch_v1.BatchServiceClient() as client:
            updated_jobs = [client.get_job(name=job.name) for job in batch_jobs]
        all_jobs_succeeded = all(
            [job.status.state == BATCH_STATE.SUCCEEDED for job in updated_jobs]
        )
        if any(
            [
                (job.status.state == BATCH_STATE.FAILED)
                or (job.status.state == BATCH_STATE.DELETION_IN_PROGRESS)
                for job in updated_jobs
            ]
        ):
            raise RuntimeError("One or more batch jobs failed")
        sleep(5)


def create_container_job(
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
    spot: bool,
    cpu_milli_per_task: int | None,
    memory_mb_per_task: int | None,
    gcs_bucket_path: str | None,
    config_file: str | None,
) -> batch_v1.Job:
    """
    This method shows how to create a sample Batch Job that will run
    a simple command inside a container on Cloud Compute instances.

    Args:
        project_id: project ID or project number of the Cloud project you want to use.
        region: name of the region you want to use to run the job. Regions that are
            available for Batch are listed on:
              https://cloud.google.com/batch/docs/get-started#locations

    Returns:
        A job object representing the job created.
    """

    add_config_to_bucket(
        bucket_name=gcs_bucket_path,
        blob_name=f"country_configs/{os.path.basename(config_file)}",
        local_file_name=config_file,
    )

    job_name = (
        "-".join(
            command.lower().replace(".", "-").replace("_", "-").replace("/", "-").split(" ")[3:4]
        )[27:]
        + f"-{datetime.now().strftime('%Y%m%d%H%M%S')}"
    )

    # Define what will be done as part of the job.
    runnable = batch_v1.Runnable()
    runnable.container = batch_v1.Runnable.Container()
    runnable.container.image_uri = f"{image}:{image_tag}"
    # NOTE: No entrypoint with micromamba docker image:
    # https://micromamba-docker.readthedocs.io/en/latest/quick_start.html#activating-a-conda-environment-for-entrypoint-commands
    runnable.container.commands = command.split(" ")
    runnable.container.options = "--ipc=host --net=host --shm-size 32G"

    # Jobs can be divided into tasks. In this case, we have only one task.
    task = batch_v1.TaskSpec()
    task.runnables = [runnable]

    # Jobs can use an existing Cloud Storage Bucket as a storage volume.
    gcs_bucket = batch_v1.GCS()
    gcs_bucket.remote_path = gcs_bucket_path
    gcs_volume = batch_v1.Volume()
    gcs_volume.gcs = gcs_bucket
    # TODO: mount volume in `ro` mode, by default it is mounting in `rw` mode.
    # NOTE: /mnt/disks/... required as mount point for volume in `rw` mode:
    # https://www.googlecloudcommunity.com/gc/Infrastructure-Compute-Storage/Seeing-new-error-mounting-GCS-bucket-on-Google-Cloud-Batch/m-p/491851
    gcs_volume.mount_path = "/mnt/disks/gcs"
    task.volumes = [gcs_volume]

    # We can specify what resources are requested by each task.
    resources = batch_v1.ComputeResource()
    resources.cpu_milli = cpu_milli_per_task or MACHINE_TYPE_TO_CPU[machine_type]
    resources.memory_mib = memory_mb_per_task or MACHINE_TYPE_TO_RAM[machine_type]
    task.compute_resource = resources

    task.max_retry_count = max_retries
    task.max_run_duration = max_duration

    # Tasks are grouped inside a job using TaskGroups.
    # Currently, it's possible to have only one task group.
    group = batch_v1.TaskGroup()
    group.task_count = task_count
    group.parallelism = parallelism
    group.task_spec = task

    # Policies are used to define on what kind of virtual machines the tasks will run on.
    policy = batch_v1.AllocationPolicy.InstancePolicy()
    policy.machine_type = machine_type
    policy.provisioning_model = 2 if spot else 0  # SPOT
    instances = batch_v1.AllocationPolicy.InstancePolicyOrTemplate()
    instances.policy = policy
    allocation_policy = batch_v1.AllocationPolicy()
    allocation_policy.instances = [instances]

    job = batch_v1.Job()
    job.task_groups = [group]
    job.allocation_policy = allocation_policy
    job.labels = {"env": "testing", "type": "container"}
    # We use Cloud Logging as it's an out of the box available option
    job.logs_policy = batch_v1.LogsPolicy()
    job.logs_policy.destination = batch_v1.LogsPolicy.Destination.CLOUD_LOGGING

    create_request = batch_v1.CreateJobRequest()
    create_request.job = job
    create_request.job_id = job_name
    # The job's parent is the region in which the job will run
    create_request.parent = f"projects/{project_id}/locations/{region}"
    
    logger.info(f"Starting job '{job_name}' with command '{command }'")
    with batch_v1.BatchServiceClient() as client:
        return client.create_job(create_request)


def add_config_to_bucket(
    bucket_name: str, blob_name: str, local_file_name: str
) -> None:
    """
    Uses the gcs_file_utils module to upload a config file from the configs
    folder into the gcs bucket.

    Args:
        project_id (str): The ID of the project.
        gcs_bucket_path (str): The path to the Google Cloud Storage bucket.
    """
    try:
        upload_file_to_bucket(
            bucket_name=bucket_name,
            blob_name=blob_name,
            local_file_name=local_file_name,
            content_type="application/x-yaml",
        )
    except Exception as e:
        print(f"An error occurred: {e}")
