#!/usr/bin/env bash

# Define the ISO codes you want to loop over
all_iso_codes=("GL")

# Set some other batch job parameters
bucket="feo-pypsa-staging"

# Loop over all iso codes
for iso in "${all_iso_codes[@]}"; do 
    # Run the python submit-job command for each country
    python submit_job \
    --command "snakemake --cores 8 $bucket/networks/$iso/elec_s_1_ec_lcopt_1H.nc --configfile /mnt/disks/gcs/$bucket/country_configs/config.$iso.yaml" \
    --image "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image" --image-tag "latest" \
    --gcs-bucket-path $bucket \
    --configfile ./ClimateTRACE/configs/config.$iso.yaml \
    --project-id "tz-feo-staging" \
    --region "europe-west2" \
    --machine-type "n1-highmem-16" \
    --disk-size-gb 128 \
    --max-duration "36000s" \
    &
done

wait
