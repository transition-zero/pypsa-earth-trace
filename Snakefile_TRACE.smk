'''

Snakefile_TRACE.smk

    These are a set of custom rules we apply for the ClimateTRACE workflow.

'''

rule ClimateTRACE:
    params:
        solving=config["solving"],
        augmented_line_connection=config["augmented_line_connection"],
    input:
        #network="feo-pypsa-staging/networks/" + RDIR + "elec_s_10_ec_lv1.00_1H.nc"
        network=GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc")
    output:
        #network="feo-pypsa-staging/results/" + RDIR + "trace-output/" + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
        network=GS.remote(BUCKET + "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc")
    log:
        solver=normpath(
            GS.remote(
                BUCKET
                + "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
            )
        ),
        python=GS.remote(
            BUCKET 
            + "logs/" 
            + RDIR 
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log"
        ),
    script:
        "scripts-custom/solve.py"