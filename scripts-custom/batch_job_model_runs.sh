# Define the ISO codes you want to loop over
iso_codes=("FJ" "FO" "TO" "GD" "KY" "GL" "TC" "KI" "VG" "MQ" "CY" "CK" "MS" "GI" "NE" "TL" "NR" "BM" "FK" "VI" "VC" "PM" "VU" "LR")
# Loop over all iso codes
for iso in "${iso_codes[@]}"
do
    # Run the python submit-job command for each country
    python submit_job --command "snakemake --cores 4 feo-pypsa-staging/networks/$iso/elec_s_1_ec_lcopt_1H.nc --configfile /mnt/disks/gcs/country_configs/config.$iso.yaml" --image "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image" --image-tag "latest" --gcs-bucket-path "feo-pypsa-staging" --config-file ./country_configs/config.$iso.yaml --project-id "tz-feo-staging" --region "europe-west2" --machine-type "n1-standard-8"
done