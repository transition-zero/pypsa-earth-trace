# -*- coding: utf-8 -*-
import asyncio
import os

import pandas as pd
from google.cloud import storage


def get_gcs_client():
    return storage.Client()


def extract_iso2_codes_from_configs(config_directory: str) -> list:
    """
    Extracts ISO2 codes from configuration files in the specified directory.

    Args:
        config_directory (str): The directory path where the configuration files are located.

    Returns:
        list: A list of ISO2 codes extracted from the configuration files.
    """
    iso2_codes = []
    for filename in os.listdir(config_directory):
        if filename.startswith("config.") and filename.endswith(".yaml"):
            iso2 = filename.split(".")[1]
            iso2_codes.append(iso2)
    return iso2_codes


async def get_demand_profile_sum(
    bucket_name: str,
    iso2: str,
    client: storage.Client,
):
    """
    Calculates the sum of demand profiles for a given ISO2 code.

    Args:
        bucket_name (str): The name of the bucket.
        iso2 (str): The ISO2 code for the country.
        client (storage.Client): The storage client.

    Returns:
        float: The total sum of demand profiles.

    Raises:
        Exception: If there is an error processing the demand profiles, returns 500.
    """

    loop = asyncio.get_event_loop()

    def read_csv():
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(f"resources/{iso2}/demand_profiles.csv")
        return pd.read_csv(blob.open("r"))

    try:
        df = await loop.run_in_executor(None, read_csv)
        df = df.drop(columns=["time"])
        total_sum = df.sum().sum()
        return total_sum
    except Exception as e:
        print(f"Error processing {iso2}: {e}")
        return 500


async def create_summary_dataframe(
    bucket_name: str, isos: list, client: storage.Client
) -> pd.DataFrame:
    """
    Creates a summary dataframe with country ISO codes and their corresponding
    total demand sums.

    Args:
        bucket_name (str): The name of the bucket.
        isos (list): A list of country ISO codes.
        client (storage.Client): The storage client.

    Returns:
        pd.DataFrame: A dataframe containing the country ISO codes and their
        total demand sums.
    """
    data = {"Country_ISO2": [], "Total_Sum": []}
    tasks = []

    for iso in isos:
        print(f"Getting demand profile for {iso}")
        tasks.append(get_demand_profile_sum(bucket_name, iso, client))

    results = await asyncio.gather(*tasks)

    for iso, total_sum in zip(isos, results):
        data["Country_ISO2"].append(iso)
        data["Total_Sum"].append(total_sum)
        print(f"calculated total sum for {iso}")

    summary_df = pd.DataFrame(data)
    return summary_df


def write_to_excel(df: pd.DataFrame, file_path: str) -> None:
    df.to_excel(file_path, index=False)


async def main(config_directory: str, bucket_name: str, output_excel_path: str) -> None:
    """
    Main function to generate a summary dataframe and write it to an Excel
    file.

    Args:
        config_directory (str): The directory containing the configuration files.
        bucket_name (str): The name of the bucket.
        output_excel_path (str): The path to the output Excel file.

    Returns:
        None
    """
    client = get_gcs_client()
    iso2_codes = extract_iso2_codes_from_configs(config_directory)
    summary_df = await create_summary_dataframe(
        bucket_name=bucket_name, isos=iso2_codes, client=client
    )
    write_to_excel(summary_df, output_excel_path)


if __name__ == "__main__":
    CONFIG_DIRECTORY = "./ClimateTrace/configs/"
    BUCKET_NAME = "feo-pypsa-staging"
    OUTPUT_EXCEL_PATH = "demand_profile_summary.xlsx"

    asyncio.run(main(CONFIG_DIRECTORY, BUCKET_NAME, OUTPUT_EXCEL_PATH))
