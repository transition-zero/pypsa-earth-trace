# -*- coding: utf-8 -*-
import logging
import os

from google.api_core.exceptions import NotFound
from google.cloud import storage

logging.basicConfig(level=logging.INFO)


def download_file_from_bucket(bucket_name: str, blob_name: str, local_file_name: str):
    """
    Downloads a file from a Google Cloud Storage bucket and loads it.

    Args:
        bucket_name (str): The name of the bucket.
        blob_name (str): The name of the blob (file) in the bucket.
        local_file_name (str): The name of the local file where the blob will be downloaded.

    Returns:
        bool: True if the file was downloaded successfully, False otherwise.
    """
    if os.path.exists(local_file_name):
        logging.info(f"File '{local_file_name}' already exists.")
        return False

    try:
        storage_client = storage.Client()
        bucket = storage_client.get_bucket(bucket_name)
        blob = bucket.blob(blob_name)
        blob.download_to_filename(local_file_name)
        return True
    except NotFound as e:
        logging.error(f"An error occurred: {e}")
        return False


def list_blobs(bucket_name: str):
    """
    Lists all the blobs in the bucket.

    Args:
        bucket_name (str): The name of the bucket.
    """

    storage_client = storage.Client()

    # Note: Client.list_blobs requires at least package version 1.17.0.
    blobs = storage_client.list_blobs(bucket_name)

    # Note: The call returns a response only when the iterator is consumed.
    for blob in blobs:
        logging.info(blob.name)


def upload_file_to_bucket(bucket_name: str, blob_name: str, local_file_name: str):
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
        blob.upload_from_filename(local_file_name)
        print(f"File {local_file_name} uploaded to {blob_name}.")
    except Exception as e:
        print(f"An error occurred: {e}")


def check_file_exists(bucket_name: str, blob_name: str):
    """
    Check if a file exists in a Google Cloud Storage bucket.

    Args:
        bucket_name (str): The name of the GCS bucket.
        blob_name (str): The name of the blob (file) in the GCS bucket.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    storage_client = storage.Client()

    bucket = storage_client.get_bucket(bucket_name)
    blob = bucket.blob(blob_name)
    # Check if the file exists in the bucket

    return blob.exists()


def copy_folder(
    folder: str,
    source_bucket_name: str = "feo-pypsa-staging",
    destination_bucket_name: str = "dispatch-data-from-tz",
):
    """
    Copy a folder from one bucket to another.

    Parameters
    ----------
    source_bucket_name : str
        The name of the source GCS bucket.
    destination_bucket_name : str
        The name of the destination GCS bucket.
    folder : str
        The path to the folder in the source bucket.
    """
    # Create a storage client
    storage_client = storage.Client()

    # Get the source and destination buckets
    source_bucket = storage_client.get_bucket(source_bucket_name)
    destination_bucket = storage_client.get_bucket(destination_bucket_name)

    # Add a trailing slash if not present
    if not folder.endswith("/"):
        folder += "/"

    # Iterate over all blobs in the folder in the source bucket
    for blob in source_bucket.list_blobs(prefix=folder):
        # Create a new blob in the destination bucket
        new_blob = destination_bucket.blob(blob.name)
        # Copy the blob from the source bucket to the new blob
        new_blob.rewrite(blob)
