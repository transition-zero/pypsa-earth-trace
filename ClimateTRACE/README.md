# PyPSA-Earth-ClimateTRACE

## Submit job on Google Cloud Batch

1. See the CLI help by running:

   ```bash
   python ClimateTRACE/scripts/submit_job --help
   ```

2. Example job submission command:

   ```bash
   python ClimateTRACE/scripts/submit_job --command 'snakemake --cores 4 feo-pypsa-staging/networks/SN/elec_s_1_ec_lcopt_1H.nc --configfile /mnt/disks/gcs/feo-pypsa-staging/country_configs/config.SN.yaml' --configfile './ClimateTRACE/configs/config.SN.yaml' --gcs-bucket-path 'feo-pypsa-staging' --image 'europe-west2-docker.pkg.dev/tz-feo-staging/feo-pypsa/pypsa-earth-image' --image-tag 'latest' --project-id 'tz-feo-staging' --region 'europe-west2' --machine-type 'n1-standard-8' --disk-size-gb 128
   ```

   This can edited as required for execution of your workflow.

## Build container on Google Cloud Build

1. Use gcloud cli to asynchronously submit a Cloud Build job for the container defined in `./infra/Dockerfile`

   ```bash
   gcloud builds submit --async --region=europe-west2 --config=./infra/cloudbuild.yaml
   ```
