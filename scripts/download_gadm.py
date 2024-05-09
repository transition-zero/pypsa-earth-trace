# -*- coding: utf-8 -*-
import os
import shutil

import requests
from _helpers import (
    configure_logging,
    create_logger,
    sets_path_to_root,
    three_2_two_digits_country,
    two_2_three_digits_country,
    two_digits_2_name_country,
)
from gcs_file_utils import (
    check_file_exists,
    download_file_from_bucket,
    upload_file_to_bucket,
)

sets_path_to_root("pypsa-earth-trace")


logger = create_logger(__name__)


def download_GADM(country_code, filepath, update=False):
    """
    Download gpkg file from GADM for a given country code.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files
    update : bool
        Update = true, forces re-download of files

    Returns
    -------
    gpkg file per country
    """
    GADM_filename = get_GADM_filename(country_code)
    GADM_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/{GADM_filename}.gpkg"
    bucket_name = "feo-pypsa-staging"
    GADM_inputfile_gpkg = os.path.join(
        bucket_name,
        "data",
        "gadm",
        GADM_filename,
        GADM_filename + ".gpkg",
    )

    # Replace with your bucket name
    download_file_from_url(GADM_url, GADM_inputfile_gpkg)
    # Check if the file exists in the bucket
    # if check_file_exists(bucket_name, GADM_inputfile_gpkg):
    #     if not update:
    #         # If the file exists in the bucket and we don't want to update it, download it to the local file path
    #         download_file_from_bucket(bucket_name, GADM_inputfile_gpkg, filepath)
    #     else:
    #         # If the file exists in the bucket and we want to update it, download it from the GADM_url
    #         download_file_from_url(GADM_url, GADM_inputfile_gpkg, out_logging)
    # else:
    #     # If the file does not exist in the bucket, download it from the GADM_url
    #     download_file_from_url(GADM_url, GADM_inputfile_gpkg, out_logging)

    #     # Upload the file to the bucket
    #     upload_file_to_bucket(bucket_name, GADM_inputfile_gpkg, GADM_inputfile_gpkg)

    return GADM_inputfile_gpkg, GADM_filename


def download_file_from_url(url, filepath):
    # if out_logging:
    #     logger.warning(
    #         f"Stage 5 of 5: {os.path.basename(filepath)} does not exist, downloading to {filepath}"
    #     )
    #  create data/osm directory
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    try:
        r = requests.get(url, stream=True, timeout=300)
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
        raise Exception(
            f"GADM server is down at {url}. Data needed for building shapes can't be extracted.\n\r"
        )
    except Exception as exception:
        raise Exception(
            f"An error happened when trying to load GADM data by {url}.\n\r"
            + str(exception)
            + "\n\r"
        )
    else:
        with open(filepath, "wb") as f:
            shutil.copyfileobj(r.raw, f)


def get_GADM_filename(country_code):
    """
    Function to get the GADM filename given the country code.
    """
    special_codes_GADM = {
        "XK": "XKO",  # kosovo
        "CP": "XCL",  # clipperton island
        "SX": "MAF",  # sint maartin
        "TF": "ATF",  # french southern territories
        "AX": "ALA",  # aland
        "IO": "IOT",  # british indian ocean territory
        "CC": "CCK",  # cocos island
        "NF": "NFK",  # norfolk
        "PN": "PCN",  # pitcairn islands
        "JE": "JEY",  # jersey
        "XS": "XSP",  # spratly
        "GG": "GGY",  # guernsey
        "UM": "UMI",  # united states minor outlying islands
        "SJ": "SJM",  # svalbard
        "CX": "CXR",  # Christmas island
    }

    if country_code in special_codes_GADM:
        return f"gadm41_{special_codes_GADM[country_code]}"
    else:
        return f"gadm41_{two_2_three_digits_country(country_code)}"


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("download_gadm")
        sets_path_to_root("pypsa-earth-trace")
    configure_logging(snakemake)
    update = snakemake.params.build_shape_options["update_file"]

    countries_list = snakemake.params.countries
    out = snakemake.output
    for country_code in countries_list:
        file_gpkg, name_file = download_GADM(country_code, out, update)
