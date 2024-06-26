# -*- coding: utf-8 -*-
"""
Description.
"""

import os
import sys

import pandas as pd
import pypsa
import trace_constraints

# imports from pypsa-earth repo
sys.path.append("scripts/")
from _helpers import configure_logging
from solve_network import prepare_network


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to
    ``pypsa.linopf.network_lopf``.

    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config

    # constrain annual generation by technology
    if "constr" in opts:
        trace_constraints.constrain_annual_generation(
            n, iso=config["countries"][0], techs=["coal", "gas", "nuclear"], year=2019
        )


def solve_network(n, config, opts="", **kwargs):
    solver_options = config["solving"]["solver"].copy()
    solver_name = solver_options.pop("name")

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    from pypsa.linopf import network_lopf

    network_lopf(
        n,
        solver_name=solver_name,
        solver_options=solver_options,
        extra_functionality=extra_functionality,
        **kwargs,
    )

    return n


if __name__ == "__main__":

    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake(
            "solve_network_TRACE",
            simpl="",
            clusters="10",
            ll="v1.25",
            opts="3H",
        )

    configure_logging(snakemake)

    tmpdir = snakemake.params.solving.get("tmpdir")
    # if tmpdir is not None:
    #     Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split("-")
    solve_opts = snakemake.params.solving["options"]

    n = pypsa.Network(snakemake.input.network)
    if snakemake.params.augmented_line_connection.get("add_to_snakefile"):
        n.lines.loc[n.lines.index.str.contains("new"), "s_nom_min"] = (
            snakemake.params.augmented_line_connection.get("min_expansion")
        )

    n = prepare_network(n, solve_opts)
    #n.export_to_netcdf('notebooks/network.nc')

    n = solve_network(
        n,
        config=snakemake.config,
        opts=opts,
        solver_dir=tmpdir,
        solver_logfile=snakemake.log.solver,
    )

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
