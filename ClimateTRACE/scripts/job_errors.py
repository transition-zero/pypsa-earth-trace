import argparse
import re
from datetime import timedelta

from google.cloud import batch, logging
from ruamel.yaml import YAML

PROJECT_ID = "tz-feo-staging"
REGION = "europe-west2"

parser = argparse.ArgumentParser()
parser.add_argument("job_name", type=str)
args = parser.parse_args()

batch_client = batch.BatchServiceClient()
logging_client = logging.Client()

# From the batch job name we can get the job uid and the job creation time
job = batch_client.get_job(name=f"projects/{PROJECT_ID}/locations/{REGION}/jobs/{args.job_name}")

# Get all the log entries that contain the error message "Error in rule"
# TODO: test other error messages such as "Missing input files"
entries = logging_client.list_entries(
    filter_=(
        f'logName="projects/{PROJECT_ID}/logs/batch_task_logs" '
        f'labels.job_uid="{job.uid}" '
        f'timestamp>="{job.create_time.isoformat()}" '
        "severity>=DEFAULT "
        '"Error in rule" '
    ),
)

# Extract the timestamp of each "Error in rule" log entry and create timestamp filters for each one
# spanning an extra 100ms to capture all log entries for the "Error in rule" message
timestamps = []
queries = []
for entry in entries:
    start = entry.timestamp
    end = start + timedelta(milliseconds=100)
    timestamps.append(
        f'(timestamp>="{start.isoformat(timespec="milliseconds")}" '
        f'AND timestamp<"{end.isoformat(timespec="milliseconds")}") '
    )
    # Also save the query string which can be used to filter the failed task in logs explorer
    labels = " ".join([f'labels.{k}="{v}"' for k, v in entry.labels.items()])
    queries.append(
        f'logName="{entry.log_name}" {labels} timestamp>="{job.create_time.isoformat()}" severity>=DEFAULT'
    )

# Do a second pass to get the full "Error in rule" failure message for each failed task
entries = logging_client.list_entries(
    filter_=(
        f'logName="projects/{PROJECT_ID}/logs/batch_task_logs" '
        f'labels.job_uid="{job.uid}" '
        f'{"OR ".join(timestamps)}'
        "severity>=DEFAULT "
    )
)

# Dump all the log entries to a string and extract all the "Error in rule" messages
log_dump = "\n".join([entry.payload for entry in entries])
error_pat = re.compile(
    r"Error in rule (?P<rule>\w+):\n\s+jobid: .*\n(?:\s+input: .*\n)?(?:\s+output: .*\n)?(?:\s+log: .*\n)?"
)
errors = error_pat.finditer(log_dump)

# Parse each "Error in rule" messages for the rule name and the run name and save these alongside
# the query string
log_info = []
for error, query in zip(errors, queries, strict=True):
    rule = error.group("rule")
    run = re.search(r"[A-Z]{2}/[0-9]{4}", error.group(0)).group(0)
    log_info.append({"run": run, "rule": rule, "query": query})

# Save the error log summary info to a yaml file
with open(f"{args.job_name}.log", "w") as f:
    yaml = YAML()
    yaml.width = 40  # NOTE: fix width to make query string easier to select for copy-pasting
    yaml.dump(log_info, f)
