#!/usr/bin/env bash

lspat="gs://feo-pypsa-staging/*/??/**"
# list all objects and remove entries which are directories only
gcloud storage ls $lspat | sed '/\/$/d' > feo-pypsa-staging.txt
# remove entries which are in 2023 and 2019 run folders
sed '/\/2023\//d; /\/2019\//d' feo-pypsa-staging.txt > sources.txt
# replace run folder name with just iso code to iso code and 2013
sed -E 's/\/([A-Z]{2})\//\/\1\/2013\//g' sources.txt > destinations.txt

while IFS= read -r source && IFS= read -r dest <&3; do
    echo "gcloud storage mv $source $dest"
    # gcloud storage mv $source $dest
done < sources.txt 3< destinations.txt
