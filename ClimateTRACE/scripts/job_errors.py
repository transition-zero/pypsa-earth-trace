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
        f'timestamp>="{job.create_time.isoformat()}" AND timestamp<="{job.update_time.isoformat()}" '
        "severity>=DEFAULT "
        '"Error in rule" '
        # '"Job stats:" OR "Error in rule" '
    ),
)

# Extract the timestamp of each "Error in rule" log entry and create timestamp filters for each one
# spanning an extra 100ms to capture all log entries for the "Error in rule" message
timestamps = []
error_queries = []
for entry in entries:
    start = entry.timestamp
    end = start + timedelta(milliseconds=100)
    timestamps.append(
        f'(timestamp>="{start.isoformat(timespec="milliseconds")}" '
        f'AND timestamp<"{end.isoformat(timespec="milliseconds")}") '
    )
    # Save the query string for jobs which errored so that they can be easily found in logs explorer
    if "Error in rule" in entry.payload:
        labels = " ".join([f'labels.{k}="{v}"' for k, v in entry.labels.items()])
        error_queries.append(
            f'logName="{entry.log_name}" {labels} '
            f'timestamp>="{job.create_time.isoformat()}" AND timestamp<="{job.update_time.isoformat()}" '
            "severity>=DEFAULT"
        )

# Do a second pass to get the full "Error in rule" failure message for each failed task
# TODO: refactor this into batches iterating over chunks of the timestamps list, there is a 20,000
# word limit on the filter
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
rule_error_pat = re.compile(
    r"Error in rule (?P<rule>\w+):\n\s+jobid: .*\n(?:\s+input: .*\n)?(?:\s+output: .*\n)?(?:\s+log: .*\n)?"
)
job_stats_pat = re.compile(r"job\s+count\n[-\s]+\n(?:\w+\s+\d+\n)+total\s+\d+\n")

# job_stats = job_stats_pat.finditer(log_dump)
rule_errors = rule_error_pat.finditer(log_dump)

# Parse each "Error in rule" messages for the rule name and the run name and save these alongside
# the query string
log_errors = []
for error, query in zip(rule_errors, error_queries, strict=True):
    rule = error.group("rule")
    run = re.search(r"[A-Z]{2}/[0-9]{4}", error.group(0)).group(0)
    log_errors.append({"run": run, "rule": rule, "query": query})

# Save the error log summary info to a yaml file
with open(f"{args.job_name}.log", "w") as f:
    yaml = YAML()
    yaml.width = 40  # NOTE: fix width to make query string easier to select for copy-pasting
    yaml.dump(log_errors, f)
