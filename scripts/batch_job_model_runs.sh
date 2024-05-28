# Define the ISO codes you want to loop over
iso_codes=("UZ" "NG" "MM" "GH" "MD" "GT" "CG" "TJ" "RE" "MU" "GP" "IS" "YE" "SD" "SV" "NI" )
# Loop over all iso codes
for iso in "${iso_codes[@]}"
do
    # Run the python submit-job command for each country
    python submit-job --command "snakemake --cores 4 feo-pypsa-staging/networks/$iso/elec_s_1_ec_lcopt_1H.nc --configfile /mnt/disks/gcs/country_configs/config.$iso.yaml" --image "europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image" --image-tag "latest" --gcs-bucket-path "feo-pypsa-staging" --config-file ./country_configs/config.$iso.yaml --project-id "tz-feo-staging" --region "europe-west2" --machine-type "n1-standard-8"
done