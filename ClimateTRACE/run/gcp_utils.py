from time import sleep
from typing import Sequence

import google.api_core.exceptions
from google.cloud import batch, compute, storage

BATCH_STATE = batch.JobStatus.State


def upload_file_to_bucket(
    bucket_name: str, blob_name: str, local_file_name: str, content_type: str = None
):
    """
    Uploads a file to a Google Cloud Storage bucket.

    Args:
        bucket_name (str): The name of the bucket.
        blob_name (str): The name of the blob (file) in the bucket.
        local_file_name (str): The name of the local file to upload.
    """

    try:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket(bucket_name)
        blob = bucket.blob(blob_name)
        blob.upload_from_filename(local_file_name, content_type=content_type)
        print(f"File {local_file_name} uploaded to {blob_name}.")
    except Exception as e:
        print(f"An error occurred: {e}")


def machine_compute_resource(*, machine_type: str, region: str, project_id: str) -> tuple[int, int]:
    """Get the CPU and memory of a machine type in a Google Cloud Compute Engine region.

    Args:
        machine_type (str): Google Cloud Compute Engine machine type.
        region (str): Compute Engine region.
        project_id (str): Google Cloud project ID.

    Returns:
        tuple[int, int]: Tuple of CPU in milli and memory in MiB.
    """
    region_zones_client = compute.RegionZonesClient()
    region_zones_resource_list = region_zones_client.list(region=region, project=project_id)
    zones = [zone.name for zone in region_zones_resource_list]

    machine_type_client = compute.MachineTypesClient()
    for zone in sorted(zones):
        try:
            machine_type_resource = machine_type_client.get(
                machine_type=machine_type,
                zone=zone,
                project=project_id,
            )
            print(f"machine type {machine_type} found in {zone}")
            break
        except google.api_core.exceptions.NotFound:
            machine_type_resource = None
            continue

    if machine_type_resource is None:
        raise ValueError(f"Machine type {machine_type} not found in any zone in region {region}")

    print(
        f"machine type {machine_type} has {machine_type_resource.guest_cpus} vCPUs "
        f"and {machine_type_resource.memory_mb} MB memory"
    )
    cpu_milli = int(machine_type_resource.guest_cpus * 1000)
    memory_mib = int(machine_type_resource.memory_mb)  # * 10**6 / 2**20 (MB to MiB conversion)
    return cpu_milli, memory_mib


def wait_for_jobs_to_succeed(batch_jobs: Sequence[batch.Job]):
    """
    Block until all Batch jobs have succeeded. Raise exception if any job
    fails.

    Args:
        batch_jobs (Sequence[batch_v1.Job]): The Batch jobs to wait for.
    """
    all_jobs_succeeded = False
    while all_jobs_succeeded is False:
        with batch.BatchServiceClient() as client:
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
    image_uri: str,
    commands: list[str],
    parallelism: int,
    machine_type: str,
    spot: bool,
    disk_size_gb: int | None = None,
    gcs_bucket_path: str | None = None,
    task_environments: list[dict[str, str]] | None = None,
    job_id: str | None = None,
) -> batch.Job:
    """
    This method shows how to create a sample Batch Job that will run a simple
    command inside a container on Cloud Compute instances.

    Args:
        project_id: project ID or project number of the Cloud project you want to use.
        region: name of the region you want to use to run the job. Regions that are
            available for Batch are listed on:
              https://cloud.google.com/batch/docs/get-started#locations

    Returns:
        A job object representing the job created.
    """
    # Jobs can be divided into tasks. In this case, we have only one task.
    task = batch.TaskSpec()

    # Jobs can use an existing Cloud Storage Bucket as a storage volume.
    if gcs_bucket_path is not None:
        gcs_bucket = batch.GCS()
        gcs_bucket.remote_path = gcs_bucket_path
        gcs_volume = batch.Volume()
        gcs_volume.gcs = gcs_bucket
        gcs_volume.mount_path = f"/mnt/disks/gcs/{gcs_bucket_path}"
        task.volumes = [gcs_volume]

    # Define what will be done as part of the job.
    runnable = batch.Runnable()
    runnable.container = batch.Runnable.Container()
    runnable.container.image_uri = image_uri
    runnable.container.commands = commands
    # TODO: use secret_variables for the license file instead of mounting it
    runnable.container.options = f"--env GRB_LICENSE_FILE={gcs_volume.mount_path}/gurobi.lic"
    task.runnables = [runnable]

    # We can specify what resources are requested by each task.
    resources = batch.ComputeResource()
    resources.cpu_milli, resources.memory_mib = machine_compute_resource(
        machine_type=machine_type, region=region, project_id=project_id
    )
    task.compute_resource = resources

    # Tasks are grouped inside a job using TaskGroups.
    # Currently, it's possible to have only one task group.
    group = batch.TaskGroup()
    group.parallelism = parallelism
    group.task_count_per_node = 1
    group.task_spec = task
    if task_environments is not None:
        group.task_environments = [batch.Environment(variables=env) for env in task_environments]

    # Policies are used to define on what kind of virtual machines the tasks will run on.
    policy = batch.AllocationPolicy.InstancePolicy()
    policy.machine_type = machine_type
    policy.provisioning_model = 2 if spot else 0  # SPOT
    if disk_size_gb is not None:
        boot_disk = batch.AllocationPolicy.Disk()
        boot_disk.size_gb = disk_size_gb
        policy.boot_disk = boot_disk
    instances = batch.AllocationPolicy.InstancePolicyOrTemplate()
    instances.policy = policy
    allocation_policy = batch.AllocationPolicy()
    allocation_policy.instances = [instances]

    job = batch.Job()
    job.task_groups = [group]
    job.allocation_policy = allocation_policy
    # We use Cloud Logging as it's an out of the box available option
    job.logs_policy = batch.LogsPolicy()
    job.logs_policy.destination = batch.LogsPolicy.Destination.CLOUD_LOGGING

    create_request = batch.CreateJobRequest()
    create_request.job = job
    if job_id is not None:
        create_request.job_id = job_id
    # The job's parent is the region in which the job will run
    create_request.parent = f"projects/{project_id}/locations/{region}"

    with batch.BatchServiceClient() as client:
        return client.create_job(create_request)
