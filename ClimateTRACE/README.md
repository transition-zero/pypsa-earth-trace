# PyPSA-Earth-ClimateTRACE

## Run Snakemake workflow on Google Cloud Batch

1. See the CLI help by running:

   ```bash
   python ClimateTRACE/run --help
   ```

2. Example job submission command:

   ```shell
   python ClimateTRACE/run --target=solve_all_networks --iso-include SZ MM BF ZW --weather-years 2019 2023 --snakemake-extra-args="--rerun-triggers mtime --cores 1" --batch-mode=tasks --task-parallelism=2 --machine-type=n1-standard-8 --disk-size-gb=64 
   ```

   This can edited as required for execution of your workflow.

## Build container on Google Cloud Build

1. Use gcloud cli to asynchronously submit a Cloud Build job for the container defined in `infra/Dockerfile`

   ```bash
   gcloud builds submit --config=infra/cloudbuild.yaml --region=europe-west2 --async
   ```
