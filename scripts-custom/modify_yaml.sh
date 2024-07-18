#!/bin/bash

# Directory containing the config files
CONFIG_DIR="./country_configs"

# Define the arguments for the Python script
YEAR_START=2023
LOAD_SCALE=1.0
ENABLE_RULES=true

# Iterate through all config.{ISO}.yaml files in the directory
for config_file in $CONFIG_DIR/config.*.yaml; do
    # Check if the file exists
    if [[ -f "$config_file" ]]; then
        echo "Processing $config_file"

        # Run the Python script with the specified arguments
        python3 update_yaml.py --file_path "$config_file" --year_start $YEAR_START --load_scale $LOAD_SCALE --enable_rules
    fi
done
