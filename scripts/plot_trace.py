import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from typing import Union
import os
from _helpers import configure_logging, create_logger, get_ember_data
import country_converter as coco

CARRIER_COLOURS = {
    "Bioenergy": "#467850",
    "Coal": "#000000",
    "Gas": "#9ABBCA",
    "Geothermal": "#d05094",
    "Nuclear": "#5F605C",
    "Oil": "#1D565C",
    "Wind": "#1F82C0",
    "Run of River": "#79A8CA",
    "Solar": "#FFD744",
    "load": "#98df8a",
    "Lignite": "#553834",
}


def set_plot_style():
    plt.style.use(
        [
            "classic",
            "seaborn-white",
            {
                "axes.grid": False,
                "grid.linestyle": "--",
                "grid.color": "0.6",
                "hatch.color": "white",
                "patch.linewidth": 0.5,
                "font.size": 12,
                "legend.fontsize": "medium",
                "lines.linewidth": 1.5,
                "pdf.fonttype": 42,
            },
        ]
    )


def generation_data(
    n: pypsa.Network, aggregate_time: Union[str, bool]
) -> Union[pd.DataFrame, pd.Series]:
    """
    Extracts generation data from the network.

    Parameters:
    n (pypsa.Network): The network object.
    aggregate_time (Union[str, bool]): The aggregation method for time.
    Can be 'mean', 'sum', or False.

    Returns:
    pd.DataFrame: The generation data.
    """
    if aggregate_time == "sum":
        divisor = 1e6
    else:
        divisor = 1e3
    generation_data = (
        n.statistics.supply(comps=["Generator"], aggregate_time=aggregate_time)
        .droplevel(0)
        .drop(index=["load"])
        .div(divisor)
        .T.rename(
            lambda x: "Bioenergy" if x == "Biomass" else ("Gas" if "Gas" in x else x)
        )
    )

    wind_technologies = generation_data.filter(like="Wind", axis=0)
    wind_sum = wind_technologies.sum()
    generation_data = generation_data.drop(index=wind_technologies.index)
    generation_data.loc["Wind"] = wind_sum
    return generation_data


def prepare_annual_ember_data(iso3: str, FilePath: str) -> pd.DataFrame:
    """
    Prepare annual ember data for a specific country and year.

    Args:
        iso3 (str): The ISO3 country code.
        FilePath (str): The file path to the CSV file containing the data.

    Returns:
        pd.DataFrame: A DataFrame containing the filtered and sorted data.

    """
    filter = (
        f"`Country code` == '{iso3}' and "
        '`Category` == "Electricity generation" and '
        '`Subcategory` == "Fuel" and '
        "`Year` == 2023 and "
        '`Unit` == "TWh"'
    )

    df = (
        pd.read_csv(FilePath)
        .query(filter)
        .loc[:, ["Variable", "Value"]]
        .set_index("Variable")
        .sort_values(by="Value", ascending=False)
    )
    return df


def plot_annual_gen(
    generation_data: pd.Series,
    ember_data: pd.DataFrame,
    carrier_colours: dict,
) -> plt.subplots:
    """
    Plot the annual generation data for different carriers.

    Parameters:
        generation_data (pd.Series): A pandas Series containing the generation data for each carrier.
        ember_data (pd.DataFrame): A pandas DataFrame containing the ember data for each carrier.
        carrier_colours (dict): A dictionary mapping carriers to their respective colors.

    Returns:
        fig: The matplotlib figure object.

    """
    fig, ax = plt.subplots()
    bottom = 0
    bottom_new = 0
    sorted_generation_data = generation_data.sort_values(ascending=False)
    # Plot each carrier as a stacked bar
    for carrier, value in sorted_generation_data.items():
        ax.bar(
            "Modelled",
            value,
            bottom=bottom,
            color=carrier_colours.get(carrier, "#000000"),
            label=carrier,
        )
        bottom += value
    for carrier, value in ember_data.iterrows():
        ax.bar(
            "Ember",
            value,
            bottom=bottom_new,
            color=carrier_colours.get(carrier, "#000000"),
        )
        bottom_new += value
    ax.set_title("Annual Generation in TWh", pad=20)
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left", title="Carriers")
    set_plot_style()
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "plot_trace",
            simpl="",
            clusters="10",
            ll="v1.25",
            opts="1H",
        )
    configure_logging(snakemake)

    opts = snakemake.wildcards.opts.split("-")
    country = snakemake.params.country
    iso3 = coco.convert(names=country, to="ISO3")
    n = pypsa.Network(snakemake.input[0])
    annual_data = generation_data(n, aggregate_time="sum")
    ember_annual = prepare_annual_ember_data(iso3=iso3, FilePath=snakemake.input[1])
    fig = plot_annual_gen(annual_data, ember_annual, CARRIER_COLOURS)
    fig.savefig(snakemake.output[0])
