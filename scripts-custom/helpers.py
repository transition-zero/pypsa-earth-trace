# -*- coding: utf-8 -*-
import os

import country_converter as coco
import numpy as np
import pandas as pd
from google.cloud import storage


def upload_blob(bucket_name, source_file_name, destination_blob_name):
    """
    Uploads a file to the Google Cloud Storage bucket.
    """
    # Initialize a client
    storage_client = storage.Client()

    # Get the bucket
    bucket = storage_client.bucket(bucket_name)

    # Create a blob object
    blob = bucket.blob(destination_blob_name)

    # Upload the file
    blob.upload_from_filename(source_file_name)

    # print(f'File {source_file_name} uploaded to {destination_blob_name}.')


def download_blob(bucket_name, source_blob_name, destination_file_name):
    """
    Downloads a blob from the Google Cloud Storage bucket.
    """
    # Initialize a client
    storage_client = storage.Client()

    # Get the bucket
    bucket = storage_client.bucket(bucket_name)

    # Get the blob
    blob = bucket.blob(source_blob_name)

    # Download the blob to a local file
    blob.download_to_filename(destination_file_name)

    # print(f'File {source_blob_name} downloaded to {destination_file_name}.')


def get_country_networks_from_bucket(bucket_name="feo-pypsa-staging", country_iso="MX"):

    # initiate storage client
    storage_client = storage.Client()

    # get bucket
    bucket = storage_client.get_bucket(bucket_name)

    network_files = [
        f"results/{country_iso}/networks/elec_s_10_ec_lv1.00_1H-constr_trace.nc",
        f"results/{country_iso}/networks/elec_s_10_ec_lv1.00_1H_trace.nc",
    ]

    for file in network_files:
        blob = bucket.blob(file)
        if not os.path.exists(f"../feo-pypsa-staging/results/{country_iso}/"):
            os.makedirs(f"../feo-pypsa-staging/results/{country_iso}/")

        # save locally
        f = file.replace("networks/", "")
        blob.download_to_filename(f"../feo-pypsa-staging/{f}")

        print(f"Downloaded file: {file}")


def get_modelling_progress(
    bucket_name="feo-pypsa-staging",
    save=True,
) -> pd.DataFrame:
    """
    Returns a DataFrame containing information on our modelling progress.
    """
    # initiate storage client
    storage_client = storage.Client()

    # get bucket
    bucket = storage_client.get_bucket(bucket_name)

    # get all iso codes for which we have pre-solved networks

    iso_codes_networks = list(
        set(
            [
                i.name.replace("networks/", "")[0:2]
                for i in bucket.list_blobs()
                if "elec_s" in i.name and "networks/" in i.name
            ]
        )
    )

    # get all iso codes for which we have results for the unconstrained scenario
    iso_codes_results_unconstr = list(
        set(
            [
                i.name.replace("results/", "")[0:2]
                for i in bucket.list_blobs()
                if "results" in i.name and "1H" in i.name
            ]
        )
    )

    # get all iso codes for which we have results for the constrained scenario
    iso_codes_results_constr = list(
        set(
            [
                i.name.replace("results/", "")[0:2]
                for i in bucket.list_blobs()
                if "results" in i.name and "1H-constr" in i.name
            ]
        )
    )

    # make pandas dataframe
    progress_data = pd.read_csv("_TRACE_outputs/model-benchmarks.csv")[
        [
            "iso",
            "country",
            "total_pwr_emissions_historic_ember_mtco2_2019",
            "total_pwr_emissions_share_historic_ember_%_2019",
        ]
    ].dropna()

    # get all iso codes for which we have results for the unconstrained scenario
    iso_codes_results_unconstr = list(
        set(
            [
                i.name.replace("results/", "")[0:2]
                for i in bucket.list_blobs()
                if "results" in i.name and "1H" in i.name
            ]
        )
    )

    # get all iso codes for which we have results for the constrained scenario
    iso_codes_results_constr = list(
        set(
            [
                i.name.replace("results/", "")[0:2]
                for i in bucket.list_blobs()
                if "results" in i.name and "1H-constr" in i.name
            ]
        )
    )

    # make pandas dataframe
    progress_data = pd.read_csv("_TRACE_outputs/model-benchmarks.csv")[
        [
            "iso",
            "country",
            "total_pwr_emissions_historic_ember_mtco2_2019",
            "total_pwr_emissions_share_historic_ember_%_2019",
        ]
    ].dropna()

    progress_data.columns = [
        "iso",
        "country",
        "total_emissions_2019",
        "share_of_global_pwr_emissions_2019",
    ]

    # assign columns
    progress_data.loc[
        progress_data.iso.isin(iso_codes_networks), "presolve_networks"
    ] = True
    progress_data.loc[
        progress_data.iso.isin(iso_codes_results_unconstr), "solved_model_unconstr"
    ] = True
    progress_data.loc[
        progress_data.iso.isin(iso_codes_results_constr), "solved_model_annual_matching"
    ] = True

    progress_data = progress_data.fillna(False).sort_values(
        by=["presolve_networks", "solved_model_unconstr"], ascending=False
    )

    # ensure data is in the correct format
    progress_data.total_emissions_2019 = progress_data.total_emissions_2019.astype(
        float
    ).round(2)
    progress_data.share_of_global_pwr_emissions_2019 = (
        progress_data.share_of_global_pwr_emissions_2019.astype(float).round(2)
    )

    # add pypsa-eur countries
    # https://github.com/PyPSA/pypsa-eur/blob/master/config/config.default.yaml

    pypsa_eur_iso = [
        "AL",
        "AT",
        "BA",
        "BE",
        "BG",
        "CH",
        "CZ",
        "DE",
        "DK",
        "EE",
        "ES",
        "FI",
        "FR",
        "GB",
        "GR",
        "HR",
        "HU",
        "IE",
        "IT",
        "LT",
        "LU",
        "LV",
        "ME",
        "MK",
        "NL",
        "NO",
        "PL",
        "PT",
        "RO",
        "RS",
        "SE",
        "SI",
        "SK",
    ]

    # for col in [
    #     "presolve_networks",
    #     "solved_model_unconstr",
    #     "solved_model_annual_matching",
    # ]:
    #     progress_data.loc[progress_data.iso.isin(pypsa_eur_iso), col] = True

    if save:
        progress_data.to_csv("_TRACE_outputs/model-progress.csv", index=False)

        upload_blob(
            bucket_name="feo-pypsa-staging",
            source_file_name="_TRACE_outputs/model-progress.csv",
            destination_blob_name="_TRACE_output_tables/model-progress.csv",
        )

    # print progress
    print('')
    print('************************************************************************************************')
    print('SUMMARY STATISTICS:')
    print('')
    sum_covered = progress_data.loc[
        (progress_data.solved_model_unconstr == True)
        & (progress_data.solved_model_annual_matching == True)
    ].share_of_global_pwr_emissions_2019.sum()

    # share of emissions covered
    print(
        f"> We are currently covering {round(sum_covered,1)}% of global power sector emissions in 2019"
    )

    # number of countries with presolved networks
    print(
        f'> Number of countries with presolved networks: { progress_data.loc[progress_data.presolve_networks == True].shape[0] }'      
    )

    # number of countries with solved models (unconstrained)
    print(
        f'> Number of countries with solved models (unconstrained): { progress_data.loc[progress_data.solved_model_unconstr == True].shape[0] }'      
    )

    # number of countries with solved models (annual matching constraint)
    print(
        f'> Number of countries with solved models (with annual matching): { progress_data.loc[progress_data.solved_model_annual_matching == True].shape[0] }'      
    )

    print('')
    print('************************************************************************************************')
    print('')

    return progress_data


def get_all_prepared_networks_from_bucket(bucket_name="feo-pypsa-staging") -> list:
    """
    Get a list of successful runs from the bucket.
    """
    # initiate storage client
    storage_client = storage.Client()
    # get bucket
    bucket = storage_client.get_bucket(bucket_name)
    # return list
    return list(
        set(
            [
                i.name.replace("networks/", "")[0:2]
                for i in bucket.list_blobs()
                if "networks" in i.name
            ]
        )
    )


def get_ember_data(
    api_key,
    dataset="electricity-generation",
    resolution="monthly",
) -> pd.DataFrame:

    import requests

    print("Loading data from Ember API... this can take a while!")

    # Define the base URL of the API
    base_url = "https://api.ember-climate.org/"
    endpoint = f"v0/{dataset}/{resolution}"

    # Make a GET request to fetch the data with headers
    response = requests.get(f"{base_url}{endpoint}?api_key={api_key}")

    print(f"API response status: {response.status_code}")

    data = pd.DataFrame(response.json()["data"])

    # change iso codes from three letter to two letter
    if "entity_code" in data.columns:
        three_letter_iso = data["entity_code"].unique()
        two_letter_iso = coco.convert(three_letter_iso, to="ISO2")
        iso_mapping = {
            three_letter_iso[i]: two_letter_iso[i] for i in range(len(three_letter_iso))
        }
        data["entity_code"] = data["entity_code"].map(iso_mapping)

    # remove world data
    if "World" in data.entity.unique():
        data.loc[data.entity_code == "not found", "entity_code"] = np.nan

    # remap variable names
    if "series" in data.columns:

        tech_mapping = {
            "Bioenergy": "biomass",
            "Coal": "coal",
            "Gas": "gas",
            "Hydro": "hydro",
            "Nuclear": "nuclear",
            "Other Fossil": np.nan,
            "Other Renewables": np.nan,
            "Solar": "solar",
            "Wind": "wind",
        }

        # overwrite existing variable names
        data["series"] = data["series"].map(tech_mapping)

    # save
    data.to_csv(f"../data/ember-{dataset}-{resolution}.csv", index=False)

    return data.dropna(axis=0, inplace=False)


def make_tracker_sheet(
    save=False,
    base_year=2019,
):

    historical = get_historical_data(
        path_to_data="../data/ember_electricity_data.csv",
        field_to_get="Power sector emissions",
    )

    # get iso codes for dataframe
    iso_codes = historical.index.get_level_values(0).unique()

    # make dataframe
    trace_tracker = pd.DataFrame(
        {"iso": iso_codes, "country": coco.convert(names=iso_codes, to="name_short")}
    )

    # append emissions
    trace_tracker[
        f"total_pwr_emissions_historic_ember_mtco2_{base_year}"
    ] = trace_tracker["iso"].map(
        historical.xs(base_year, level=1)
        .reset_index()
        .groupby(by=["Country code"])
        .sum(numeric_only=True)
        .Value.to_dict()
    )

    # add emissions share
    trace_tracker[f"total_pwr_emissions_share_historic_ember_%_{base_year}"] = (
        trace_tracker[f"total_pwr_emissions_historic_ember_mtco2_{base_year}"]
        / trace_tracker[f"total_pwr_emissions_historic_ember_mtco2_{base_year}"].sum()
        * 100
    )

    # append coal generation
    coal_generation = (
        get_historical_data(path_to_data="../data/ember_electricity_data.csv")
        .query('Unit == "TWh"')
        .query('Variable == "coal"')
        .xs(base_year, level=1)
        .Value.to_dict()
    )

    trace_tracker[f"coal_generation_historic_ember_TWh_{base_year}"] = (
        trace_tracker["iso"].map(coal_generation).fillna(0)
    )

    # append trace estimates
    trace_estimates = (
        pd.read_csv("../data/power/electricity-generation_emissions-sources.csv")
        .dropna(subset=["activity", "activity_units"], axis=0)
        # .query(' source_type.str.contains("coal") ')
        .query(' source_type == "coal" ')
        .query(f' start_time.str.contains("{base_year}") ')
        .groupby(["iso3_country"], as_index=False)
        .sum(numeric_only=True)
    )

    trace_estimates = (
        trace_estimates.assign(
            iso=coco.convert(trace_estimates["iso3_country"], to="iso2")
        )[["iso", "activity"]]
        .set_index("iso")
        .activity.divide(1e6)  # convert to TWh
        .to_dict()
    )

    trace_tracker[f"coal_generation_historic_TRACE_TWh_{base_year}"] = (
        trace_tracker["iso"].map(trace_estimates).fillna(0)
    )

    # to csv
    if save:
        trace_tracker.to_csv("../_TRACE_outputs/model-benchmarks.csv", index=False)

    return trace_tracker
