import pandas as pd
import os
import pycountry

csv_file_path_ember = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "annual_demand_totals_ember.csv"
)

csv_file_path_trace = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "demand_profile_summary.xlsx"
)


def iso3_to_iso2(iso3):
    try:
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
trace = pd.read_excel(csv_file_path_trace).rename(columns={"Country_ISO2": "iso2"})


def calculate_scale_factor(row):
    """
    Calculate the scale factor for a row.
    """
    if row["Value_MWh"] != 0 and row["Total_Sum"] != 0:
        return row["Total_Sum"] / row["Value_MWh"]
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
    factors = produce_factors(ember, trace, year=2019)
    factors.to_csv("demand_scale_factors.csv", index=False)
