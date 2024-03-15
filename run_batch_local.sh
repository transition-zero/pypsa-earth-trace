#!/bin/bash

# ------------------------------------------------------------
#   run_batch_local.sh
#
#   This script allows you to run multiple network solves
#    with the PyPSA-earth workflow in batch on your local
#    machine. While PyPSA-earth has built-in functionality
#    to do this to some extent, we wrote this script because 
#    we want to similar workflows run many countries for TRACE.
#
#   @amanmajid (13/03/2024)
#
# ------------------------------------------------------------

echo "------------------------------------------------------------------------------------------"
echo "BEGIN"
echo ""

# Define the list of codes
#   For now we run just one iso code at a time because of restrictions on our Gurobi license
iso_codes=("MX")
number_of_clusters=("10") #"50"
opts=("3H")

# Loop through each iso code
for iso_code in "${iso_codes[@]}"; do
    # loop through each cluster number
    for cluster in "${number_of_clusters[@]}"; do
        # loop through hourly resolution
        for opt in "${opts[@]}"; do

            echo ""
            echo ">>>>"
            echo "Running ${iso_code} with ${cluster} clusters and ${hour} hourly resolution"
            echo "<<<<"
            echo ""

            # run snakecommand
            snakemake -c5 -j5 \
            feo-pypsa-staging/results/${iso_code}/trace-output/elec_s_${cluster}_ec_lv1.25_${opt}.nc \
            --configfile country_configs/config.${iso_code}.yaml \
            --snakefile Snakefile \
            # --dry-run \
            # -n \
            # -F

        done
    done
done

echo ""
echo "------------------------------------------------------------------------------------------"

#snakemake -c5 -j5 feo-pypsa-staging/networks/MX/elec_s_25_ec_lv1.25_1H.nc --configfile country_configs/config.MX.yaml 