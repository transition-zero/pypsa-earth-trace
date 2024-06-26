"""

Snakefile_TRACE.smk

    These are a set of custom rules we apply for the ClimateTRACE workflow.

"""

rule solve_all_networks_TRACE:
    input:
        GS.remote(
            expand(
                BUCKET
                + "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc",
                **config["scenario"],
            )
        ),

rule solve_network_TRACE:
    params:
        solving=config["solving"],
        augmented_line_connection=config["augmented_line_connection"],
    input:
        network=GS.remote(
            BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
        ),
    output:
        network=GS.remote(
            BUCKET
            + "results/"
            + RDIR
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc"
        ),
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
