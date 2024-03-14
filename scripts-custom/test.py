
"""
Description
"""

import os
import sys
import pypsa
import pandas as pd

sys.path.append("scripts/")
from _helpers import configure_logging

def demand_calibration():
    pass

def constrain_annual_generation():
    pass

def constrain_monthly_generation():
    pass

def constrain_monthly_generation():
    pass

def do_something(input_path,output_path):
    """
    Simple test.
    """
    print('>>>>>>>>>>>>>>>>>>>>>')
    print('RUN SCRIPT: TEST.PY')
    print('<<<<<<<<<<<<<<<<<<<<<')

    n = pypsa.Network(input_path)
    n.export_to_netcdf(output_path)
        

if __name__ == "__main__":

    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("modify_demand")

    configure_logging(snakemake)

    do_something(
        input_path=snakemake.input.input_network,
        output_path=snakemake.output.output_network,
    )