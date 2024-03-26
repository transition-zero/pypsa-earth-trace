import pypsa
import pandas as pd

import helpers

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
        n : pypsa.Network,
        iso : str,
        techs : list,
        year : int,
    ) -> None:

    # get historical data
    try:
        historic_data = pd.read_csv('data/ember_electricity_data.csv')
    except:
        ValueError('Could not find historical data at data/ember_electricity_data.csv')

    for i, tech in enumerate(techs):

        # calculate historic value
        iso = iso
        year = str(year)

        historical_value = (
            historic_data
            .query(" series == @tech ")
            .query(" entity_code == @iso ")
            .query(" date.str.contains(@year)")
            .groupby(by='entity_code')
            .sum(numeric_only=True)
            .generation_twh
            .mul(1e6)
            .round(0)
            .iloc[0]
        )

        # get unique generators
        gens_i = n.generators.query(' carrier == @tech ').index

        # define lhs
        lhs_gen = \
            linexpr(
                (1, get_var(n, "Generator", "p")[gens_i].T)
            ).sum().sum()
        
        # define rhs
        rhs = historical_value

        print(f'Constraining annual generation for {tech} to {rhs:.2f} MWh')

        define_constraints(n, lhs_gen, "==", rhs, f"annual_{tech}_gen")

def constrain_monthly_generation():
    pass