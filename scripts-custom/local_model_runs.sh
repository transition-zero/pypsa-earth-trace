#!/bin/bash

# ------------------------------------------------------------
#   local_model_runs.sh
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
    #"1H-constr"
)

# Define whether we want transmission to be expandable
# - 1.00 means that transmission is fixed
# - 1.25 means that transmission can be expanded by 25%
ll="v1.00"


# ----------------------------------------
# Option 1:
# - Loop over all iso codes

# Navigate to the directory
cd ../country_configs || exit
# Declare an empty array to store iso codes
iso_codes=()
# Iterate over each file in the directory
for i in config.*; do
    # Extract the two remaining letters
    iso=$(echo "$i" | cut -d '.' -f 2)
    # Add iso to the array
    iso_codes+=("$iso")
done

#----------------------------------------
# Option 2:
# - Run over a single iso code only

#iso_codes=("VN")
#iso_codes=('SL' 'GA' 'GF' 'MT' 'GN' 'MU' 'MR' 'MG' 'CG' 'LB' 'CF' 'SS' 'SR' 'SG' 'SN' 'ML' 'BH' 'NC' 'GY' 'GW' 'BI' 'AF' 'LS' 'TT' 'TG')

iso_codes=(
"AS"
"AG"
"AW"
"BS"
"BB"
"BZ"
"BJ"
"BT"
"BF"
"CV"
"KY"
"TD"
"KM"
"CK"
"DJ"
"DM"
"ER"
"FK"
"FO"
"PF"
"GM"
"GI"
"GD"
"GP"
"GU"
"HT"
"LR"
"MV"
"MQ"
"MS"
"NR"
"NE"
"NU"
"PG"
"KN"
"LC"
"VC"
"WS"
"ST"
"SC"
"SB"
"SO"
"TO"
"TC"
"VU"
"VG"
"VI"
"EH"
)

# change director
cd .. || exit

# ---
# Run models

for iso_code in "${iso_codes[@]}"; do
    # loop through each cluster number
    for cluster in "${number_of_clusters[@]}"; do
        # loop through hourly resolution
        for opt in "${opts[@]}"; do

            # Reduce the number of clusters for countries with smaller grids

            # Countries with a single bus
            max_cluster=1
            if [[ "$iso_code" == "SL" || \
                  "$iso_code" == "GA" || \
                  "$iso_code" == "MT" || \
                  "$iso_code" == "CG" || \
                  "$iso_code" == "CF" || \
                  "$iso_code" == "SS" || \
                  "$iso_code" == "SR" || \
                  "$iso_code" == "ML" || \
                  "$iso_code" == "BH" || \
                  "$iso_code" == "NC" || \
                  "$iso_code" == "GY" || \
                  "$iso_code" == "SB" || \
                  "$iso_code" == "SC" || \
                  "$iso_code" == "VU" || \
                  "$iso_code" == "VI" || \
                  "$iso_code" == "VG" || \
                  "$iso_code" == "MS" || \
                  "$iso_code" == "KM" || \
                  "$iso_code" == "NR" || \
                  "$iso_code" == "NU" || \
                  "$iso_code" == "BB" || \
                  "$iso_code" == "GI" || \
                  "$iso_code" == "AW" || \
                  "$iso_code" == "KY" || \
                  "$iso_code" == "GU" || \
                  "$iso_code" == "BJ" || \
                  "$iso_code" == "LC" || \
                  "$iso_code" == "FK" || \
                  "$iso_code" == "TC" || \
                  "$iso_code" == "ST" || \
                  "$iso_code" == "SO" || \
                  "$iso_code" == "VC" || \
                  "$iso_code" == "PF" || \
                  "$iso_code" == "AG" || \
                  "$iso_code" == "HT" || \
                  "$iso_code" == "BF" || \
                  "$iso_code" == "KN" || \
                  "$iso_code" == "DJ" || \
                  "$iso_code" == "BS" || \
                  "$iso_code" == "DM" || \
                  "$iso_code" == "NE" || \
                  "$iso_code" == "GM" || \
                  "$iso_code" == "AS" || \
                  "$iso_code" == "GD" || \
                  "$iso_code" == "MV" || \
                  "$iso_code" == "CV" || \
                  "$iso_code" == "FO" || \
                  "$iso_code" == "CK" || \
                  "$iso_code" == "TO" || \
                  "$iso_code" == "WS" || \
                  "$iso_code" == "GW" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 2 buses
            max_cluster=2
            if [[ "$iso_code" == "BI" || \
                  "$iso_code" == "LS" || \
                  "$iso_code" == "TD" || \
                  "$iso_code" == "PG" || \
                  "$iso_code" == "LR" || \
                  "$iso_code" == "TT" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 3 buses
            max_cluster=3
            if [[ "$iso_code" == "GF" || \
                  "$iso_code" == "MR" || \
                  "$iso_code" == "MG" || \
                  "$iso_code" == "LB" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 4 buses
            max_cluster=4
            if [[ "$iso_code" == "SN" || \
                  "$iso_code" == "AF" || \
                  "$iso_code" == "ER" || \
                  "$iso_code" == "TG" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 5 buses
            max_cluster=5
            if [[ "$iso_code" == "SG" || \
                  "$iso_code" == "BZ" || \
                  "$iso_code" == "MQ" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi
            
            # Countries with 6 buses
            max_cluster=6
            if [[ "$iso_code" == "GN" || \
                  "$iso_code" == "MU" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 7 buses
            max_cluster=7
            if [[ "$iso_code" == "EH" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 8 buses
            max_cluster=8
            if [[ "$iso_code" == "CM" || \
                  "$iso_code" == "BT" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            # Countries with 9 buses
            max_cluster=9
            if [[ "$iso_code" == "GP" ]]; then
                if (( cluster > max_cluster )); then
                    cluster=$max_cluster
                fi
            fi

            echo ""
            echo "Running ${iso_code} with ${cluster} clusters"
            echo ""

            # run snakecommand
            snakemake -c2 -j2 \
            feo-pypsa-staging/results/${iso_code}/networks/elec_s_${cluster}_ec_l${ll}_${opt}_trace.nc \
            --configfile country_configs/config.${iso_code}.yaml \
            --snakefile Snakefile #--unlock #\
            #-F
            #--unlock
            #--rerun-incomplete

        done
    done
done

# run benchmarks
# python scripts-custom/benchmarking.py

echo ""
echo "------------------------------------------------------------------------------------------"