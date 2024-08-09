# PyPSA-Earth-ClimateTRACE

## Run Snakemake workflow on Google Cloud Batch

1. See the CLI help by running:

   ```bash
   python ClimateTRACE/run --help
   ```

2. Example job submission command:

   ```bash
   python ClimateTRACE/run --target=solve_all_networks --weather-year=2023 --snakemake-extra-args='--rerun-triggers mtime --cores 1' --batch-mode=tasks --task-parallelism=2 --machine-type=n1-standard-8 --disk-size-gb=64 --iso-include SZ MM BF ZW
   ```

   This can edited as required for execution of your workflow.

## Build container on Google Cloud Build

1. Use gcloud cli to asynchronously submit a Cloud Build job for the container defined in `infra/Dockerfile`

   ```bash
   gcloud builds submit --config=infra/cloudbuild.yaml --region=europe-west2 --async
   ```
