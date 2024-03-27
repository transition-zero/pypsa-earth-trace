#!/bin/bash

# ------------------------------------------------------------
#   run_models_in_batch.sh
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

# Define number of clusters (nodes)
#  - For multiple arguments, input as a list such as ("10" "20" "30")
number_of_clusters=("10") 

# Define options
# - 1H runs the model without constraints
# - 1H-constr runs the model with annual matching
# - For multiple arguments, input as a list such as ("1H-constr" "1H")
opts=(
    "1H-constr"
    #"1H-constr" "1H"
) 

# Define whether we want transmission to be expandable
# - 1.00 means that transmission is fixed
# - 1.25 means that transmission can be expanded by 25%
ll="v1.00"

# Get all iso codes we can loop over
# - We get this from the country_configs directory

# Navigate to the directory
cd country_configs || exit
# Declare an empty array to store iso codes
iso_codes=()
# Iterate over each file in the directory
for i in config.*; do
    # Extract the two remaining letters
    iso=$(echo "$i" | cut -d '.' -f 2)
    # Add iso to the array
    iso_codes+=("$iso")
done

cd .. || exit

#iso_codes=("MX")

# ---
# Run models

for iso_code in "${iso_codes[@]}"; do
    # loop through each cluster number
    for cluster in "${number_of_clusters[@]}"; do
        # loop through hourly resolution
        for opt in "${opts[@]}"; do

            echo ""
            echo "Running ${iso_code} with ${cluster} clusters"
            echo ""

            # run snakecommand
            snakemake -c5 -j5 \
            feo-pypsa-staging/results/${iso_code}/networks/elec_s_${cluster}_ec_l${ll}_${opt}_trace.nc \
            --configfile country_configs/config.${iso_code}.yaml \
            --snakefile Snakefile \
            -F
            #--unlock

        done
    done
done

echo ""
echo "------------------------------------------------------------------------------------------"