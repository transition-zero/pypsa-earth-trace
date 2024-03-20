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
    
    if n.config['ClimateTrace']['annual_generation']:

        # get historical data
        historical = (
            helpers
            .get_historical_data(path_to_data='data/ember_electricity_data.csv')
            .query(' Unit == "TWh" ')
        )

        for i, tech in enumerate(techs):

            # get unique generators
            gens_i = n.generators.query(' carrier == @tech ').index

            # define lhs
            lhs_gen = \
                linexpr(
                    (1, get_var(n, "Generator", "p")[gens_i].T)
                ).sum().sum()

            # define rhs
            rhs = historical.loc[(iso, year)].query(" Variable == @tech ").Value.iloc[0] * 1e6

            define_constraints(n, lhs_gen, "==", rhs, f"annual_{tech}_gen")


def constrain_monthly_generation():
    pass