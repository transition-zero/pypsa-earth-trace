# -*- coding: utf-8 -*-
import glob
import os

import matplotlib.pyplot as plt
import pypsa

# Get a list of all .nc files in the directory
files = glob.glob("../results/*/networks/*.nc")
color_mapping = {
    "onwind": "#235ebc",
    "onshore wind": "#235ebc",
    "offwind": "#6895dd",
    "offwind-ac": "#6895dd",
    "offshore wind": "#6895dd",
    "offshore wind ac": "#6895dd",
    "offwind-dc": "#74c6f2",
    "offshore wind dc": "#74c6f2",
    "hydro": "#08ad97",
    "hydro+PHS": "#08ad97",
    "PHS": "#08ad97",
    "hydro reservoir": "#08ad97",
    "hydroelectricity": "#08ad97",
    "ror": "#4adbc8",
    "run of river": "#4adbc8",
    "solar": "#f9d002",
    "solar PV": "#f9d002",
    "solar thermal": "#ffef60",
    "biomass": "#0c6013",
    "solid biomass": "#06540d",
    "biogas": "#23932d",
    "waste": "#68896b",
    "geothermal": "#ba91b1",
    "OCGT": "#d35050",
    "gas": "#d35050",
    "natural gas": "#d35050",
    "CCGT": "#b20101",
    "nuclear": "#ff9000",
    "coal": "#707070",
    "lignite": "#9e5a01",
    "oil": "#262626",
    "H2": "#ea048a",
    "hydrogen storage": "#ea048a",
    "battery": "#b8ea04",
    "Electric load": "#f9d002",
    "electricity": "#f9d002",
    "lines": "#70af1d",
    "transmission lines": "#70af1d",
    "AC-AC": "#70af1d",
    "AC line": "#70af1d",
    "links": "#8a1caf",
    "HVDC links": "#8a1caf",
    "DC-DC": "#8a1caf",
    "DC link": "#8a1caf",
    "load": "#FF0000",
    # Add more mappings as needed
}

for file in files:
    # Load the network
    n = pypsa.Network(file)
    country_iso = os.path.normpath(file).split(os.sep)[1]
    # Create the plot
    data = (
        n.statistics.dispatch(
            comps=["Generator"], aggregate_time=False, groupby="carrier"
        )
        .droplevel(0)
        .T.iloc[0:48]
    )
    colors = [color_mapping.get(label, "black") for label in data.columns]
    ax = data.plot(kind="bar", stacked=True, color=colors)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.set_title(f"{country_iso} - Generation by Carrier")
    # Extract the country ISO from the file path

    # Set the title
    title = f'{country_iso}_{os.path.basename(file).split(".")[0]}'  # Use the country ISO and filename (without extension) as the title
    ax.set_title(title)
    os.makedirs(f"results/{country_iso}/plots", exist_ok=True)
    # Save the plot as a PNG
    plt.savefig(f"results/{country_iso}/plots/{title}.png", bbox_inches="tight")

    # Close the plot to free up memory
    plt.close()
