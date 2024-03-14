'''

Snakefile_TRACE.smk

    These are a set of custom rules we apply for the ClimateTRACE workflow.

'''

rule TRACE_TEST:
    input:
        network="feo-pypsa-staging/results/" + RDIR + "networks/" + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    output:
        network="feo-pypsa-staging/results/" + RDIR + "trace-output/" + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    script:
        "scripts-custom/test.py"