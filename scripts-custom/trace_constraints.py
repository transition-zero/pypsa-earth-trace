# -*- coding: utf-8 -*-
import helpers
import pandas as pd
import pypsa

# linopf
from pypsa.linopf import (
    define_constraints,
    define_variables,
    get_var,
    join_exprs,
    linexpr,
)


def demand_calibration():
    pass


def constrain_annual_generation(
    n: pypsa.Network,
    iso: str,
    techs: list,
    year: int,
) -> None:

    # get historical data
    try:
        historic_data = pd.read_csv("data/ember-electricity-generation-yearly.csv")
    except:
        ValueError(
            "Could not load historical generation data! Please make sure you run get_ember_data() from helpers.py"
        )

    for i, tech in enumerate(techs):
        # calculate historic value
        iso = iso
        year = year
        historical_value = (
            historic_data
            .query(" series == @tech ")
            .query(" entity_code == @iso ")
            .query(" date == @year ")
            .groupby(by="entity_code")
            .sum(numeric_only=True)
            .generation_twh.mul(1e6)
            .round(0)
            .iloc[0]
        )

        # get unique generators
        if tech == "gas":
            gens_i = n.generators.query(' carrier.isin(["OCGT", "CCGT"]) ').index
        else:
            gens_i = n.generators.query(" carrier == @tech ").index

        # define lhs
        lhs_gen = linexpr((1, get_var(n, "Generator", "p")[gens_i].T)).sum().sum()

        #Set constraint only if technology exist within network
        if lhs_gen != 0:
            # define rhs
            rhs = historical_value
            # print for reference
            print(f"Constraining annual generation for {tech} to {rhs:.2f} MWh")
            # define constraint
            define_constraints(n, lhs_gen, "==", rhs, f"annual_{tech}_gen")


def constrain_monthly_generation():
    pass
