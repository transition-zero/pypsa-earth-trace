#!/usr/bin/env bash

## Define the ISO code to touch output for
iso="RU"

## Name of the bucket containing files to touch
bucket="feo-pypsa-staging"

## Touch the output files in the order they would be produced by the workflow

# build_shapes
gcloud storage mv gs://$bucket/resources/$iso/shapes gs://$bucket/resources/$iso/shapes.tmp
gcloud storage mv gs://$bucket/resources/$iso/shapes.tmp gs://$bucket/resources/$iso/shapes

# download_osm_data
gcloud storage mv gs://$bucket/resources/$iso/osm/raw gs://$bucket/resources/$iso/osm/raw.tmp
gcloud storage mv gs://$bucket/resources/$iso/osm/raw.tmp gs://$bucket/resources/$iso/osm/raw

# clean_osm_data
gcloud storage mv gs://$bucket/resources/$iso/osm/clean gs://$bucket/resources/$iso/osm/clean.tmp
gcloud storage mv gs://$bucket/resources/$iso/osm/clean.tmp gs://$bucket/resources/$iso/osm/clean

# build_osm_network
gcloud storage mv gs://$bucket/resources/$iso/base_network gs://$bucket/resources/$iso/base_network.tmp
gcloud storage mv gs://$bucket/resources/$iso/base_network.tmp gs://$bucket/resources/$iso/base_network

# base_network
gcloud storage mv gs://$bucket/networks/$iso/base.nc gs://$bucket/networks/$iso/base.nc.tmp
gcloud storage mv gs://$bucket/networks/$iso/base.nc.tmp gs://$bucket/networks/$iso/base.nc

# build_bus_regions
gcloud storage mv gs://$bucket/resources/$iso/bus_regions gs://$bucket/resources/$iso/bus_regions.tmp
gcloud storage mv gs://$bucket/resources/$iso/bus_regions.tmp gs://$bucket/resources/$iso/bus_regions

# build_powerplants
gcloud storage mv gs://$bucket/resources/$iso/powerplants.csv gs://$bucket/resources/$iso/powerplants.csv.tmp
gcloud storage mv gs://$bucket/resources/$iso/powerplants.csv.tmp gs://$bucket/resources/$iso/powerplants.csv
gcloud storage mv gs://$bucket/resources/$iso/powerplants_osm2pm.csv gs://$bucket/resources/$iso/powerplants_osm2pm.csv.tmp
gcloud storage mv gs://$bucket/resources/$iso/powerplants_osm2pm.csv.tmp gs://$bucket/resources/$iso/powerplants_osm2pm.csv

# build_demand_profiles
gcloud storage mv gs://$bucket/resources/$iso/demand_profiles.csv gs://$bucket/resources/$iso/demand_profiles.csv.tmp
gcloud storage mv gs://$bucket/resources/$iso/demand_profiles.csv.tmp gs://$bucket/resources/$iso/demand_profiles.csv

# build_cutout
gcloud storage mv gs://$bucket/cutouts/$iso gs://$bucket/cutouts/$iso.tmp
gcloud storage mv gs://$bucket/cutouts/$iso.tmp gs://$bucket/cutouts/$iso

# build_renewable_profiles
gcloud storage mv gs://$bucket/resources/$iso/renewable_profiles gs://$bucket/resources/$iso/renewable_profiles.tmp
gcloud storage mv gs://$bucket/resources/$iso/renewable_profiles.tmp gs://$bucket/resources/$iso/renewable_profiles

# add_electricity
gcloud storage mv gs://$bucket/networks/$iso/elec.nc gs://$bucket/networks/$iso/elec.nc.tmp
gcloud storage mv gs://$bucket/networks/$iso/elec.nc.tmp gs://$bucket/networks/$iso/elec.nc

# logs
gcloud storage mv gs://$bucket/logs/$iso gs://$bucket/logs/$iso.tmp
gcloud storage mv gs://$bucket/logs/$iso.tmp gs://$bucket/logs/$iso

# benchmarks
gcloud storage mv gs://$bucket/benchmarks/$iso gs://$bucket/benchmarks/$iso.tmp
gcloud storage mv gs://$bucket/benchmarks/$iso.tmp gs://$bucket/benchmarks/$iso
