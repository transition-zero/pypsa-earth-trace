import pandas as pd
import os


def modify_fuels(cost_df: pd.DataFrame, fuel_df: pd.DataFrame) -> pd.DataFrame:
    """Replace standard fuel costs with those from the fuel_df based on snapshot year.
        This also needs the iso2 for that country.
    Args:
        cost_df (pd.DataFrame): _description_
        fuel_df (pd.DataFrame): _description_

    Returns:
        pd.DataFrame: _description_
    """
    return NotImplementedError


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        # TODO - add path to the fuel data file in the Snakefile input
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "retrieve_cost_data",
        )
        # Propose to load the fuel data for all the years.
        costs_df = pd.read_csv(snakemake.input[0])
        costs_df = modify_fuels(costs_df)

        costs_df.to_csv(snakemake.output[0], index=False)
