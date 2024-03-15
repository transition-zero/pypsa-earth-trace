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
        techs : list,
    ) -> None:
    
    if n.config['ClimateTrace']['annual_generation']:

        # get historical data
        historical = helpers.get_historical_data(
            path_to_data='data/ember_electricity_data.csv'
        )

        # get unique generators
        gens_i = n.generators.query("carrier in @techs").index

        # define lhs
        lhs_gen = \
            linexpr(
                (1, get_var(n, "Generator", "p")[gens_i].T)
            ).sum().sum()

        # define rhs
        #rhs = historical.loc[('MEX', 2019)].query(" Variable == 'Coal' ").Value.iloc[0] * 1e6
        rhs = 50

        define_constraints(n, lhs_gen, "==", rhs, f"{techs[0]}_gen")


def constrain_monthly_generation():
    pass