#!/usr/bin/env bash

## Name of the bucket containing files to touch
bucket="feo-pypsa-staging"

## Define scenario wildcard parameters to touch output for
simpl=""
clusters="1"
ll="v1.25"
opts="1H"

## Define run outputs to touch
if [ $# -eq 0 ]; then
  runs=($(gcloud storage ls "gs://${bucket}/resources/*/" | sed '/^$/d' | sed '/^.*:$/d' | rev | cut -d'/' -f2-3 | rev))
  runs=(${runs[@]:1})
else
  runs=($@)
fi

touch_gcloud() {
  gcloud storage mv gs://${bucket}/$1 gs://${bucket}/$1.tmp && \
  gcloud storage mv gs://${bucket}/$1.tmp gs://${bucket}/$1
}

run_touch() {
    ## Touch the output files in the order they would be produced by the workflow

    # build_shapes
    touch_gcloud resources/$1/shapes

    # build_cutout
    touch_gcloud cutouts/$1

    # download_osm_data
    touch_gcloud resources/$1/osm/raw

    # clean_osm_data
    touch_gcloud resources/$1/osm/clean

    # build_osm_network
    touch_gcloud resources/$1/base_network

    # base_network
    touch_gcloud networks/$1/base.nc

    # build_powerplants
    touch_gcloud resources/$1/powerplants.csv
    touch_gcloud resources/$1/powerplants_osm2pm.csv

    # build_bus_regions
    touch_gcloud resources/$1/bus_regions/regions_onshore.geojson
    touch_gcloud resources/$1/bus_regions/regions_offshore.geojson

    # build_demand_profiles
    touch_gcloud resources/$1/demand_profiles.csv

    # copy_defaultnatura_tiff
    touch_gcloud resources/$1/natura.tiff

    # build_renewable_profiles
    touch_gcloud resources/$1/renewable_profiles

    # retrieve_cost_data
    touch_gcloud resources/$1/costs.csv

    # modify_cost_data
    touch_gcloud resources/$1/costs_modified.csv

    # add_electricity
    touch_gcloud networks/$1/elec.nc

    # simplify_network
    touch_gcloud networks/$1/elec_s${simpl}.nc
    touch_gcloud resources/$1/bus_regions/regions_onshore_elec_s${simpl}.geojson
    touch_gcloud resources/$1/bus_regions/regions_offshore_elec_s${simpl}.geojson
    touch_gcloud resources/$1/bus_regions/busmap_elec_s${simpl}.csv
    touch_gcloud resources/$1/bus_regions/connection_costs_s${simpl}.csv

    # cluster_network
    touch_gcloud networks/$1/elec_s${simpl}_${clusters}.nc
    touch_gcloud resources/$1/bus_regions/regions_onshore_elec_s${simpl}_${clusters}.geojson
    touch_gcloud resources/$1/bus_regions/regions_offshore_elec_s${simpl}_${clusters}.geojson
    touch_gcloud resources/$1/bus_regions/busmap_elec_s${simpl}_${clusters}.csv
    touch_gcloud resources/$1/bus_regions/linemap_elec_s${simpl}_${clusters}.csv

    # add_extra_components
    touch_gcloud networks/$1/elec_s${simpl}_${clusters}_ec.nc

    # prepare_network
    touch_gcloud networks/$1/elec_s${simpl}_${clusters}_ec_l${ll}_${opts}.nc

    # solve_network
    touch_gcloud results/$1/networks/elec_s${simpl}_${clusters}_ec_l${ll}_${opts}_trace.nc

    # logs
    touch_gcloud logs/$1

    # benchmarks
    touch_gcloud benchmarks/$1
}


## Run touch script concurrently in batches
batch_size=10
mkdir -p /tmp/trace

for ((i=0; i<${#runs[@]}; i+=batch_size)); do
  batch=(${runs[@]:i:batch_size})
  echo "Touching outputs for runs: ${batch[@]}"

  for run in "${batch[@]}"; do
    log_file="/tmp/trace/$(echo $run | sed 's/\//-/g').log"
    echo "run: $run"
    run_touch $run >$log_file 2>&1 </dev/null &
  done

  ## Wait for all background processes in the batch to complete
  sleep 1
  echo "Waiting for batch to complete..."
  wait

done
