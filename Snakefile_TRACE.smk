'''

Snakefile_TRACE.smk

    These are a set of custom rules we apply for the ClimateTRACE workflow.

'''

rule TRACE_TEST:
    input:
        input_network='feo-pypsa-staging/results/MX/networks/elec_s_10_ec_lv1.25_1H.nc'
    output:
        output_network='feo-pypsa-staging/results/MX/networks/elec_s_10_ec_lv1.25_1H_TEST.nc'
    script:
        "scripts-custom/test.py"