import os
import sys

import pandas as pd
import pycountry
from pandas._libs.parsers import STR_NA_VALUES

try:
    STR_NA_VALUES.remove("NA")
    keep_default_na = False
    na_values = STR_NA_VALUES
except KeyError:
    keep_default_na = True
    na_values = None

csv_file_path_ember = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "../trace_data/annual_demand_totals_ember.csv",
)

csv_file_path_trace = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "../trace_data/demand_profile_summary.xlsx",
)


def iso3_to_iso2(iso3):
    try:
        if iso3 == "XKX":
            return "XK"
        return pycountry.countries.get(alpha_3=iso3).alpha_2
    except AttributeError:
        return None


def convert_twh_to_mwh(df):
    df["Value_MWh"] = df["Value"] * 1000000
    return df


ember = (
    pd.read_csv(csv_file_path_ember)
    .assign(Country_code=lambda df: df["Country code"].apply(iso3_to_iso2))
    .drop(
        columns=[
            "Country code",
            "Area",
            "Area type",
            "Category",
            "Subcategory",
            "Variable",
        ]
    )
    .rename(columns={"Country_code": "iso2"})
    .pipe(convert_twh_to_mwh)
)
trace = pd.read_excel(
    csv_file_path_trace,
    keep_default_na=keep_default_na,
    na_values=na_values,
).rename(columns={"Country_ISO2": "iso2"})


def calculate_scale_factor(row):
    """
    Calculate the scale factor for a row.
    """
    if row["Value_MWh"] != 0 and row["Total_Sum"] != 0:
        return row["Value_MWh"] / row["Total_Sum"]
    else:
        return None


def add_scale_factor(df):
    """
    Add scale factor to the DataFrame.
    """
    df["scale_factor"] = df.apply(calculate_scale_factor, axis=1)
    return df


def produce_factors(ember: pd.DataFrame, trace: pd.DataFrame, year: int):
    """
    Produces demand scale factors for each country by dividing the total demand in the Ember dataset by the total demand in the TRACE dataset.
    """

    ember = ember.query("Year == @year").drop(columns=["Year"])
    merged = ember.merge(trace, on="iso2", suffixes=("_ember", "_trace")).pipe(
        add_scale_factor
    )

    return merged


if __name__ == "__main__":
    year = sys.argv[1]
    factors = produce_factors(ember, trace, year=int(year))
    factors.to_csv(
        f"ClimateTRACE/trace_data/demand_scale_factors_{year}.csv", index=False
    )
