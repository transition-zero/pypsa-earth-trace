# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Creates electric demand profile csv.

Relevant Settings
-----------------

.. code:: yaml

    load:
        scale:
        ssp:
        weather_year:
        prediction_year:
        region_load:

Inputs
------

- ``networks/base.nc``: confer :ref:`base`, a base PyPSA Network
- ``resources/bus_regions/regions_onshore.geojson``: confer :mod:`build_bus_regions`
- ``load_data_paths``: paths to load profiles, e.g. hourly country load profiles produced by GEGIS
- ``resources/shapes/gadm_shapes.geojson``: confer :ref:`shapes`, file containing the gadm shapes

Outputs
-------

- ``resources/demand_profiles.csv``: the content of the file is the electric demand profile associated to each bus. The file has the snapshots as rows and the buses of the network as columns.

Description
-----------

The rule :mod:`build_demand` creates load demand profiles in correspondence of the buses of the network.
It creates the load paths for GEGIS outputs by combining the input parameters of the countries, weather year, prediction year, and SSP scenario.
Then with a function that takes in the PyPSA network "base.nc", region and gadm shape data, the countries of interest, a scale factor, and the snapshots,
it returns a csv file called "demand_profiles.csv", that allocates the load to the buses of the network according to GDP and population.
"""
import logging
import os
from itertools import product

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import scipy.sparse as sparse
import xarray as xr
from _helpers import configure_logging, create_logger, getContinent, update_p_nom_max
from shapely.prepared import prep
from shapely.validation import make_valid

logger = create_logger(__name__)


def normed(s):
    return s / s.sum()


def get_load_paths_gegis(ssp_parentfolder, config):
    """
    Create load paths for GEGIS outputs.

    The paths are created automatically according to included country,
    weather year, prediction year and ssp scenario

    Example
    -------
    ["/data/ssp2-2.6/2030/era5_2013/Africa.nc", "/data/ssp2-2.6/2030/era5_2013/Africa.nc"]
    """
    countries = config.get("countries")
    region_load = getContinent(countries)
    weather_year = config.get("load_options")["weather_year"]
    prediction_year = config.get("load_options")["prediction_year"]
    ssp = config.get("load_options")["ssp"]

    load_paths = []
    for continent in region_load:
        load_path = os.path.join(
            ssp_parentfolder,
            str(ssp),
            str(prediction_year),
            "era5_" + str(weather_year),
            str(continent) + ".nc",
        )
        load_paths.append(load_path)

    return load_paths


def prepare_plexos_demands(path: str, countries_iso2: dict) -> pd.DataFrame:
    """
    Prepare Plexos demands by reading a CSV file, transforming the data, and returning a DataFrame.

    Args:
        countries_iso2 (dict): A dictionary mapping region codes to region names.

    Returns:
        pd.DataFrame: A DataFrame containing the prepared Plexos demands data.

    """
    df = pd.read_csv(path)
    df["Datetime"] = df["Datetime"].apply(lambda x: x if ":" in x else f"{x} 00:00")
    df["Datetime"] = pd.to_datetime(
        df["Datetime"], format="%d/%m/%Y %H:%M", errors="coerce"
    )
    df_melted = df.melt(
        id_vars=["Datetime"],
        var_name="region_code",
        value_name="Electricity demand",
    )
    df_melted["region_name"] = df_melted["region_code"].map(countries_iso2)
    df_melted.set_index("Datetime", inplace=True)
    return df_melted


def shapes_to_shapes(orig, dest):
    """
    Adopted from vresutils.transfer.Shapes2Shapes()
    """
    orig_prepped = list(map(prep, orig))
    transfer = sparse.lil_matrix((len(dest), len(orig)), dtype=float)

    for i, j in product(range(len(dest)), range(len(orig))):
        if orig_prepped[j].intersects(dest[i]):
            area = orig[j].intersection(dest[i]).area
            transfer[i, j] = area / dest[i].area

    return transfer


def build_demand_profiles(
    n,
    load_paths,
    regions,
    admin_shapes,
    countries,
    scale,
    start_date,
    end_date,
    out_path,
    plexos_demand_path: str,
):
    """
    Create csv file of electric demand time series.

    Parameters
    ----------
    n : pypsa network
    load_paths: paths of the load files
    regions : .geojson
        Contains bus_id of low voltage substations and
        bus region shapes (voronoi cells)
    admin_shapes : .geojson
        contains subregional gdp, population and shape data
    countries : list
        List of countries that is config input
    scale : float
        The scale factor is multiplied with the load (1.3 = 30% more load)
    start_date: parameter
        The start_date is the first hour of the first day of the snapshots
    end_date: parameter
        The end_date is the last hour of the last day of the snapshots

    Returns
    -------
    demand_profiles.csv : csv file containing the electric demand time series
    """
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    regions = gpd.read_file(regions).set_index("name").reindex(substation_lv_i)
    load_paths = load_paths
    # Merge load .nc files: https://stackoverflow.com/questions/47226429/join-merge-multiple-netcdf-files-using-xarray
    gegis_load = xr.open_mfdataset(load_paths, combine="nested")
    gegis_load = gegis_load.to_dataframe().reset_index().set_index("time")
    # filter load for analysed countries
    gegis_load = gegis_load.loc[gegis_load.region_code.isin(countries)]
    countries_iso2 = {
        "UG": "Uganda",
        "AF": "Afghanistan",
        "BI": "Burundi",
        "PG": "Papua New Guinea",
        "LA": "Laos",
        "XK": "Kosovo",
        "GY": "Guyana",
        "BT": "Bhutan",
        "CV": "Cape Verde",
        "GF": "French Guiana",
        "GU": "Guam",
        "DM": "Dominica",
    }
    if any(country in countries_iso2.keys() for country in countries):
        gegis_load = prepare_plexos_demands(plexos_demand_path, countries_iso2)
        gegis_load = gegis_load.loc[gegis_load.region_code.isin(countries)]

    logger.info(f"Load data scaled with scaling factor {scale}.")
    gegis_load["Electricity demand"] *= scale
    shapes = gpd.read_file(admin_shapes).set_index("GADM_ID")
    shapes["geometry"] = shapes["geometry"].apply(lambda x: make_valid(x))

    def upsample(cntry, group):
        """
        Distributes load in country according to population and gdp.
        """
        l = gegis_load.loc[gegis_load.region_code == cntry]["Electricity demand"]
        if len(group) == 1:
            return pd.DataFrame({group.index[0]: l})
        else:
            shapes_cntry = shapes.loc[shapes.country == cntry]
            transfer = shapes_to_shapes(group, shapes_cntry.geometry).T.tocsr()
            gdp_n = pd.Series(
                transfer.dot(shapes_cntry["gdp"].fillna(1.0).values), index=group.index
            )
            pop_n = pd.Series(
                transfer.dot(shapes_cntry["pop"].fillna(1.0).values), index=group.index
            )

            # relative factors 0.6 and 0.4 have been determined from a linear
            # regression on the country to EU continent load data
            # (refer to vresutils.load._upsampling_weights)
            # TODO: require adjustment for Africa
            factors = normed(0.6 * normed(gdp_n) + 0.4 * normed(pop_n))
            return pd.DataFrame(
                factors.values * l.values[:, np.newaxis],
                index=l.index,
                columns=factors.index,
            )

    demand_profiles = pd.concat(
        [
            upsample(cntry, group)
            for cntry, group in regions.geometry.groupby(regions.country)
        ],
        axis=1,
    )

    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date) - pd.Timedelta(hours=1)
    # demand_profiles = demand_profiles.loc[start_date:end_date]
    # TODO: generalise to snapshots that have a start and end date
    # which do not align with 1st January to 31st December.
    demand_profiles.index = n.snapshots
    demand_profiles.to_csv(out_path, header=True)

    logger.info(f"Demand_profiles csv file created for the corresponding snapshots.")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand_profiles")
        sets_path_to_root("pypsa-earth-trace")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.base_network)

    # Snakemake imports:
    regions = snakemake.input.regions
    load_paths = snakemake.input["load"]
    countries = snakemake.params.countries
    admin_shapes = snakemake.input.gadm_shapes
    scale = float(snakemake.params.load_options["scale"])
    start_date = snakemake.params.snapshots["start"]
    end_date = snakemake.params.snapshots["end"]
    out_path = snakemake.output[0]
    plexos_demand_path = snakemake.input.plexos_demand_path
    build_demand_profiles(
        n,
        load_paths,
        regions,
        admin_shapes,
        countries,
        scale,
        start_date,
        end_date,
        out_path,
        plexos_demand_path,
    )
