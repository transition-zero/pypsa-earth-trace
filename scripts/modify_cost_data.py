import os

import pandas as pd
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def modify_fuels(cost_df: pd.DataFrame, fuel_df: pd.DataFrame, year: int) -> pd.DataFrame:
    """Replace standard fuel costs with those from the fuel_df based on snapshot year.
        This also needs the iso2 for that country.
    Args:
        cost_df (pd.DataFrame): Standard cost data csv coming from technology data url.
        fuel_df (pd.DataFrame):fuel cost database created for the TRACE project.
    Returns:
        pd.DataFrame: Modified cost data with fuel costs replaced.
    """
    fuels = fuel_df.query(f"year == {year}")
    if fuels.empty:
        raise ValueError(f"No fuel data available for the year {year}")
    for technology in fuels.technology:
        cost_df.loc[
            (cost_df["technology"] == technology) & (cost_df["parameter"] == "fuel"),
            "value",
        ] = fuels.loc[(fuels["technology"] == technology), "value"].values

    return cost_df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "modify_cost_data",
        )
    configure_logging(snakemake)
    year = int(snakemake.output[0].split("_")[-1].split(".")[0])
    fuel_df = pd.read_csv(snakemake.input.fuel_database)
    costs_df = pd.read_csv(snakemake.input.cost_data)
    costs_df = modify_fuels(cost_df=costs_df, fuel_df=fuel_df, year=year)
    costs_df.to_csv(snakemake.output[0], index=False)
