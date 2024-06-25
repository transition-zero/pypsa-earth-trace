# Define the ISO codes you want to loop over
iso_codes=("GL")

# Set some other command parameters
machine_type="n1-highmem-16"
disk_size_gb=128
cores=8
bucket="feo-pypsa-staging"
max_duration="36000s"

# Loop over all iso codes
for iso in "${iso_codes[@]}"
do
    # Run the python submit-job command for each country
    python submit_job --command "snakemake --cores $cores $bucket/networks/$iso/elec_s_1_ec_lcopt_1H.nc --configfile /mnt/disks/gcs/$bucket/country_configs/config.$iso.yaml" --image "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image" --image-tag "latest" --gcs-bucket-path $bucket --config-file ./country_configs/config.$iso.yaml --project-id "tz-feo-staging" --region "europe-west2" --machine-type $machine_type --disk-size-gb $disk_size_gb --max-duration $max_duration
done
