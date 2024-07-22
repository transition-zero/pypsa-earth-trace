#!/usr/bin/env bash

# Define the ISO codes to run
all_iso_codes=("RU")

# Define network solve parameters
all_clusters=("3")
all_opts=("1H-constr")
ll="v1.00"
simpl=""

# Set some other batch job parameters
bucket="feo-pypsa-staging"

for iso in "${all_iso_codes[@]}"; do
    # loop through each cluster number
    for clusters in "${all_clusters[@]}"; do
        # loop through hourly resolution
        for opts in "${all_opts[@]}"; do
            # submit batch job
            python submit_job \
            --command "snakemake --cores 1 $bucket/results/$iso/networks/elec_s${simpl}_${clusters}_ec_l${ll}_${opts}_trace.nc --configfile /mnt/disks/gcs/$bucket/country_configs/config.$iso.yaml" \
            --image "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image" --image-tag "latest" \
            --gcs-bucket-path $bucket \
            --configfile ./ClimateTRACE/configs/config.$iso.yaml \
            --project-id "tz-feo-staging" \
            --region "europe-west2" \
            --machine-type "n1-standard-8" \
            --disk-size-gb 64 \
            --max-duration "36000s" \
            &
        done
    done
done

wait
