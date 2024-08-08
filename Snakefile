# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir
from shutil import copyfile, move

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

from _helpers import create_country_list, get_last_commit_message
from build_demand_profiles import get_load_paths_gegis
from retrieve_databundle_light import datafiles_retrivedatabundle
from pathlib import Path

HTTP = HTTPRemoteProvider()
GS = GSRemoteProvider(keep_local=True)

if "config" not in globals() or not config:  # skip when used as sub-workflow
    if not exists("config.yaml"):
        copyfile("config.tutorial.yaml", "config.yaml")

    configfile: "config.yaml"


configfile: "configs/bundle_config.yaml"


config.update({"git_commit": get_last_commit_message(".")})

# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])

# create a list of iteration steps, required to solve the experimental design
# each value is used as wildcard input e.g. solution_{unc}
config["scenario"]["unc"] = [
    f"m{i}" for i in range(config["monte_carlo"]["options"]["samples"])
]

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""
BUCKET = config["remote"]["gcs_bucket_path"]

load_data_paths = get_load_paths_gegis("data", config)

if config["enable"].get("retrieve_cost_data", True):
    COSTS_ = GS.remote(BUCKET + "resources/" + RDIR + "costs.csv")
else:
    COSTS_ = "data/costs.csv"

if config["enable"].get("modify_cost_data", False):
    COSTS = GS.remote(BUCKET + "resources/" + RDIR + f"costs_modified.csv")
else:
    COSTS = COSTS_

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+(m|flex|min)?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    unc="[-+a-zA-Z0-9\.]*",


if config["custom_rules"] is not []:
    for rule in config["custom_rules"]:

        include: rule


rule clean:
    run:
        try:
            shell("snakemake -j 1 solve_all_networks --delete-all-output")
        except:
            shell("snakemake -j 1 solve_all_networks_monte --delete-all-output")
            pass
        shell("snakemake -j 1 run_all_scenarios --delete-all-output")


rule run_tests:
    output:
        touch("tests.done"),
    run:
        import os

        shell("snakemake --cores all build_test_configs")
        directory = "test/tmp"  # assign directory
        for filename in os.scandir(directory):  # iterate over files in that directory
            if filename.is_file():
                print(
                    f"Running test: config name '{filename.name}'' and path name '{filename.path}'"
                )
                if "custom" in filename.name:
                    shell("mkdir -p configs/scenarios")
                    shell("cp {filename.path} configs/scenarios/config.custom.yaml")
                    shell("snakemake --cores 1 run_all_scenarios --forceall")
                if "monte" in filename.name:
                    shell("cp {filename.path} config.yaml")
                    shell("snakemake --cores all solve_all_networks_monte --forceall")
                else:
                    shell("cp {filename.path} config.yaml")
                    shell("snakemake --cores all solve_all_networks --forceall")
        print("Tests are successful.")


rule solve_all_networks:
    input:
        GS.remote(
            expand(
                BUCKET
                + "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc",
                **config["scenario"],
            )
        ),


rule plot_all_p_nom:
    input:
        expand(
            "results/"
            + RDIR
            + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_p_nom.{ext}",
            **config["scenario"],
            ext=["png", "pdf"],
        ),


rule make_all_summaries:
    input:
        expand(
            "results/"
            + RDIR
            + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
            **config["scenario"],
            country=["all"] + config["countries"],
        ),


rule plot_all_summaries:
    input:
        expand(
            "results/"
            + RDIR
            + "plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
            summary=["energy", "costs"],
            **config["scenario"],
            country=["all"] + config["countries"],
            ext=["png", "pdf"],
        ),


if config["enable"].get("retrieve_databundle", True):

    rule retrieve_databundle_light:
        params:
            countries=config["countries"],
            tutorial=config["tutorial"],
            hydrobasins_level=config["renewable"]["hydro"]["hydrobasins_level"],
        output:  #expand(directory('{file}') if isdir('{file}') else '{file}', file=datafiles)
            expand("{file}", file=datafiles_retrivedatabundle(config)),
            directory("data/landcover"),
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "retrieve_databundle.log"),
        benchmark:
            GS.remote(BUCKET + "benchmarks/" + RDIR + "retrieve_databundle_light")
        script:
            "scripts/retrieve_databundle_light.py"


if config["enable"].get("download_osm_data", True):

    rule download_osm_data:
        params:
            countries=config["countries"],
        output:
            cables=GS.remote(
                BUCKET + "resources/" + RDIR + "osm/raw/all_raw_cables.geojson"
            ),
            generators=GS.remote(
                BUCKET + "resources/" + RDIR + "osm/raw/all_raw_generators.geojson"
            ),
            generators_csv=GS.remote(
                BUCKET + "resources/" + RDIR + "osm/raw/all_raw_generators.csv"
            ),
            lines=GS.remote(
                BUCKET + "resources/" + RDIR + "osm/raw/all_raw_lines.geojson"
            ),
            substations=GS.remote(
                BUCKET + "resources/" + RDIR + "osm/raw/all_raw_substations.geojson"
            ),
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "download_osm_data.log"),
        benchmark:
            GS.remote(BUCKET + "benchmarks/" + RDIR + "download_osm_data")
        script:
            "scripts/download_osm_data.py"


rule clean_osm_data:
    params:
        crs=config["crs"],
        clean_osm_data_options=config["clean_osm_data_options"],
    input:
        cables=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/raw/all_raw_cables.geojson"
        ),
        generators=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/raw/all_raw_generators.geojson"
        ),
        lines=GS.remote(BUCKET + "resources/" + RDIR + "osm/raw/all_raw_lines.geojson"),
        substations=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/raw/all_raw_substations.geojson"
        ),
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
        offshore_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
        ),
        africa_shape=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/africa_shape.geojson"
        ),
    output:
        generators=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_generators.geojson"
        ),
        generators_csv=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_generators.csv"
        ),
        lines=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_lines.geojson"
        ),
        substations=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_substations.geojson"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "clean_osm_data.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "clean_osm_data")
    script:
        "scripts/clean_osm_data.py"


if config["build_osm_network"].get("use_kdtree", False):
    build_osm_network_script = "scripts/build_osm_network_kdtree.py"
else:
    build_osm_network_script = "scripts/build_osm_network.py"


rule build_osm_network:
    params:
        build_osm_network=config.get("build_osm_network", {}),
        countries=config["countries"],
        crs=config["crs"],
    input:
        generators=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_generators.geojson"
        ),
        lines=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_lines.geojson"
        ),
        substations=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_substations.geojson"
        ),
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
    output:
        lines=GS.remote(
            BUCKET + "resources/" + RDIR + "base_network/all_lines_build_network.csv"
        ),
        converters=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "base_network/all_converters_build_network.csv"
        ),
        transformers=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "base_network/all_transformers_build_network.csv"
        ),
        substations=GS.remote(
            BUCKET + "resources/" + RDIR + "base_network/all_buses_build_network.csv"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_osm_network.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "build_osm_network")
    script:
        build_osm_network_script


rule build_shapes:
    params:
        build_shape_options=config["build_shape_options"],
        crs=config["crs"],
        countries=config["countries"],
    input:
        # naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
        # eez='data/bundle/eez/World_EEZ_v8_2014.shp',
        # nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
        # nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
        # nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
        eez="data/eez/eez_v11.gpkg",
    output:
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
        offshore_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
        ),
        africa_shape=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/africa_shape.geojson"
        ),
        gadm_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_shapes.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "build_shapes")
    threads: 1
    resources:
        mem_mb=3096,
    script:
        "scripts/build_shapes.py"


rule base_network:
    params:
        voltages=config["electricity"]["voltages"],
        transformers=config["transformers"],
        snapshots=config["snapshots"],
        links=config["links"],
        lines=config["lines"],
        hvdc_as_lines=config["electricity"]["hvdc_as_lines"],
        countries=config["countries"],
        base_network=config["base_network"],
    input:
        osm_buses=GS.remote(
            BUCKET + "resources/" + RDIR + "base_network/all_buses_build_network.csv"
        ),
        osm_lines=GS.remote(
            BUCKET + "resources/" + RDIR + "base_network/all_lines_build_network.csv"
        ),
        osm_converters=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "base_network/all_converters_build_network.csv"
        ),
        osm_transformers=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "base_network/all_transformers_build_network.csv"
        ),
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
        offshore_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
        ),
    output:
        GS.remote(BUCKET + "networks/" + RDIR + "base.nc"),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "base_network.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "base_network")
    threads: 1
    resources:
        mem_mb=500,
    script:
        "scripts/base_network.py"


rule build_bus_regions:
    params:
        alternative_clustering=config["cluster_options"]["alternative_clustering"],
        area_crs=config["crs"]["area_crs"],
        countries=config["countries"],
    input:
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
        offshore_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
        ),
        base_network=GS.remote(BUCKET + "networks/" + RDIR + "base.nc"),
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
        ),
    output:
        regions_onshore=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/regions_onshore.geojson"
        ),
        regions_offshore=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/regions_offshore.geojson"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_bus_regions.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "build_bus_regions")
    threads: 1
    resources:
        mem_mb=1000,
    script:
        "scripts/build_bus_regions.py"


def terminate_if_cutout_exists(config=config):
    """
    Check if any of the requested cutout files exist.
    If that's the case, terminate execution to avoid data loss.
    """
    config_cutouts = [
        d_value["cutout"] for tc, d_value in config["renewable"].items()
    ] + list(config["atlite"]["cutouts"].keys())

    for ct in set(config_cutouts):
        cutout_fl = GS.remote(BUCKET + "cutouts/" + CDIR + ct + ".nc")
        if os.path.exists(cutout_fl):
            raise Exception(
                "An option `build_cutout` is enabled, while a cutout file '"
                + cutout_fl
                + "' still exists and risks to be overwritten. If this is an intended behavior, please move or delete this file and re-run the rule. Otherwise, just disable the `build_cutout` rule in the config file."
            )


if config["enable"].get("build_cutout", False):
    terminate_if_cutout_exists(config)

    rule build_cutout:
        params:
            snapshots=config["snapshots"],
            cutouts=config["atlite"]["cutouts"],
        input:
            onshore_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
            ),
            offshore_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
            ),
        output:
            GS.remote(BUCKET + "cutouts/" + CDIR + "{cutout}.nc"),
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "build_cutout/{cutout}.log"),
        benchmark:
            GS.remote(BUCKET + "benchmarks/" + RDIR + "build_cutout_{cutout}")
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        script:
            "scripts/build_cutout.py"


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        params:
            area_crs=config["crs"]["area_crs"],
        input:
            shapefiles_land="data/landcover",
            cutouts=GS.remote(
                expand(BUCKET + "cutouts/" + CDIR + "{cutouts}.nc", **config["atlite"])
            ),
        output:
            GS.remote(BUCKET + "resources/" + RDIR + "natura.tiff"),
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "build_natura_raster.log"),
        benchmark:
            GS.remote(BUCKET + "benchmarks/" + RDIR + "build_natura_raster")
        script:
            "scripts/build_natura_raster.py"


if not config["enable"].get("build_natura_raster", False):

    rule copy_defaultnatura_tiff:
        input:
            "data/natura/natura.tiff",
        output:
            GS.remote(BUCKET + "resources/" + RDIR + "natura.tiff"),
        run:
            import shutil

            shutil.copyfile(input[0], output[0])


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['version']}/outputs/costs_{config['costs']['year']}.csv",
                keep_local=True,
            ),
        output:
            COSTS_,
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "retrieve_cost_data.log"),
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


if config["enable"].get("modify_cost_data", False):

    rule modify_cost_data:
        params:
            year=int(config["snapshots"]["start"].split("-")[0]),
        input:
            cost_data=COSTS_,
            fuel_database="data/fuel_prices_db.csv",
        output:
            COSTS,
        log:
            GS.remote(BUCKET + "logs/" + RDIR + "modify_cost_data.log"),
        script:
            "scripts/modify_cost_data.py"


rule build_demand_profiles:
    params:
        snapshots=config["snapshots"],
        load_options=config["load_options"],
        countries=config["countries"],
    input:
        base_network=GS.remote(BUCKET + "networks/" + RDIR + "base.nc"),
        regions=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/regions_onshore.geojson"
        ),
        load=load_data_paths,
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
        ),
    output:
        GS.remote(BUCKET + "resources/" + RDIR + "demand_profiles.csv"),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_demand_profiles.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "build_demand_profiles")
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/build_demand_profiles.py"


rule build_renewable_profiles:
    params:
        crs=config["crs"],
        renewable=config["renewable"],
        countries=config["countries"],
        alternative_clustering=config["cluster_options"]["alternative_clustering"],
    input:
        natura=GS.remote(BUCKET + "resources/" + RDIR + "natura.tiff"),
        copernicus="data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="data/gebco/GEBCO_2021_TID.nc",
        country_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
        ),
        offshore_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/offshore_shapes.geojson"
        ),
        hydro_capacities="data/hydro_capacities.csv",
        eia_hydro_generation="data/eia_hydro_annual_generation.csv",
        powerplants=GS.remote(BUCKET + "resources/" + RDIR + "powerplants.csv"),
        regions=lambda w: (
            GS.remote(
                BUCKET + "resources/" + RDIR + "bus_regions/regions_onshore.geojson"
            )
            if w.technology in ("onwind", "solar", "hydro")
            else GS.remote(
                BUCKET + "resources/" + RDIR + "bus_regions/regions_offshore.geojson"
            )
        ),
        cutout=lambda w: GS.remote(
            BUCKET
            + "cutouts/"
            + CDIR
            + config["renewable"][w.technology]["cutout"]
            + ".nc"
        ),
    output:
        profile=GS.remote(
            BUCKET + "resources/" + RDIR + "renewable_profiles/profile_{technology}.nc"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_renewable_profile_{technology}.log"),
    benchmark:
        GS.remote(
            BUCKET + "benchmarks/" + RDIR + "build_renewable_profiles_{technology}"
        )
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 5000,
    script:
        "scripts/build_renewable_profiles.py"


rule build_powerplants:
    params:
        geo_crs=config["crs"]["geo_crs"],
        countries=config["countries"],
        gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
        alternative_clustering=config["cluster_options"]["alternative_clustering"],
        powerplants_filter=config["electricity"]["powerplants_filter"],
    input:
        base_network=GS.remote(BUCKET + "networks/" + RDIR + "base.nc"),
        pm_config=GS.remote(BUCKET + "configs/powerplantmatching_config.yaml"),
        custom_powerplants="data/custom_powerplants.csv",
        osm_powerplants=GS.remote(
            BUCKET + "resources/" + RDIR + "osm/clean/all_clean_generators.csv"
        ),
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
        ),
    output:
        powerplants=GS.remote(BUCKET + "resources/" + RDIR + "powerplants.csv"),
        powerplants_osm2pm=GS.remote(
            BUCKET + "resources/" + RDIR + "powerplants_osm2pm.csv"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "build_powerplants.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "build_powerplants")
    threads: 1
    resources:
        mem_mb=500,
    script:
        "scripts/build_powerplants.py"


rule add_electricity:
    params:
        countries=config["countries"],
        costs=config["costs"],
        conventional=config.get("conventional", {}),
        electricity=config["electricity"],
        alternative_clustering=config["cluster_options"]["alternative_clustering"],
        renewable=config["renewable"],
        length_factor=config["lines"]["length_factor"],
    input:
        **{
            f"profile_{tech}": GS.remote(
                BUCKET + "resources/" + RDIR + f"renewable_profiles/profile_{tech}.nc"
            )
            for tech in config["renewable"]
            if tech in config["electricity"]["renewable_carriers"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network=GS.remote(BUCKET + "networks/" + RDIR + "base.nc"),
        tech_costs=COSTS,
        powerplants=GS.remote(BUCKET + "resources/" + RDIR + "powerplants.csv"),
        #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
        #using this line instead of the following will test updated gadm shapes for MA.
        #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
        #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
        gadm_shapes=GS.remote(
            BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
        ),
        hydro_capacities="data/hydro_capacities.csv",
        demand_profiles=GS.remote(BUCKET + "resources/" + RDIR + "demand_profiles.csv"),
    output:
        GS.remote(BUCKET + "networks/" + RDIR + "elec.nc"),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "add_electricity.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "add_electricity")
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_electricity.py"


rule simplify_network:
    params:
        renewable=config["renewable"],
        geo_crs=config["crs"]["geo_crs"],
        cluster_options=config["cluster_options"],
        countries=config["countries"],
        build_shape_options=config["build_shape_options"],
        electricity=config["electricity"],
        costs=config["costs"],
        config_lines=config["lines"],
        config_links=config["links"],
        focus_weights=config.get("focus_weights", None),
    input:
        network=GS.remote(BUCKET + "networks/" + RDIR + "elec.nc"),
        tech_costs=COSTS,
        regions_onshore=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/regions_onshore.geojson"
        ),
        regions_offshore=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/regions_offshore.geojson"
        ),
    output:
        network=GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}.nc"),
        regions_onshore=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "bus_regions/regions_onshore_elec_s{simpl}.geojson"
        ),
        regions_offshore=GS.remote(
            BUCKET
            + "resources/"
            + RDIR
            + "bus_regions/regions_offshore_elec_s{simpl}.geojson"
        ),
        busmap=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv"
        ),
        connection_costs=GS.remote(
            BUCKET + "resources/" + RDIR + "bus_regions/connection_costs_s{simpl}.csv"
        ),
    log:
        GS.remote(BUCKET + "logs/" + RDIR + "simplify_network/elec_s{simpl}.log"),
    benchmark:
        GS.remote(BUCKET + "benchmarks/" + RDIR + "simplify_network/elec_s{simpl}")
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/simplify_network.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == True:

    rule cluster_network:
        params:
            build_shape_options=config["build_shape_options"],
            electricity=config["electricity"],
            costs=config["costs"],
            length_factor=config["lines"]["length_factor"],
            renewable=config["renewable"],
            geo_crs=config["crs"]["geo_crs"],
            countries=config["countries"],
            cluster_options=config["cluster_options"],
            focus_weights=config.get("focus_weights", None),
            #custom_busmap=config["enable"].get("custom_busmap", False)
        input:
            network=GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}.nc"),
            country_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
            ),
            regions_onshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_onshore_elec_s{simpl}.geojson"
            ),
            regions_offshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_offshore_elec_s{simpl}.geojson"
            ),
            #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
            ),
            # busmap=ancient('resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network=GS.remote(
                BUCKET
                + "networks/"
                + RDIR
                + "elec_s{simpl}_{clusters}_pre_augmentation.nc"
            ),
            regions_onshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
            busmap=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv"
            ),
            linemap=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/linemap_elec_s{simpl}_{clusters}.csv"
            ),
        log:
            GS.remote(
                BUCKET
                + "logs/"
                + RDIR
                + "cluster_network/elec_s{simpl}_{clusters}.log"
            ),
        benchmark:
            GS.remote(
                BUCKET
                + "benchmarks/"
                + RDIR
                + "cluster_network/elec_s{simpl}_{clusters}"
            )
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/cluster_network.py"

    rule augmented_line_connections:
        params:
            lines=config["lines"],
            augmented_line_connection=config["augmented_line_connection"],
            hvdc_as_lines=config["electricity"]["hvdc_as_lines"],
            electricity=config["electricity"],
            costs=config["costs"],
        input:
            tech_costs=COSTS,
            network=GS.remote(
                BUCKET
                + "networks/"
                + RDIR
                + "elec_s{simpl}_{clusters}_pre_augmentation.nc"
            ),
            regions_onshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
        output:
            network=GS.remote(
                BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}.nc"
            ),
        log:
            GS.remote(
                BUCKET
                + "logs/"
                + RDIR
                + "augmented_line_connections/elec_s{simpl}_{clusters}.log"
            ),
        benchmark:
            GS.remote(
                BUCKET
                + "benchmarks/"
                + RDIR
                + "augmented_line_connections/elec_s{simpl}_{clusters}"
            )
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/augmented_line_connections.py"


if config["augmented_line_connection"].get("add_to_snakefile", False) == False:

    rule cluster_network:
        params:
            build_shape_options=config["build_shape_options"],
            electricity=config["electricity"],
            costs=config["costs"],
            length_factor=config["lines"]["length_factor"],
            renewable=config["renewable"],
            geo_crs=config["crs"]["geo_crs"],
            countries=config["countries"],
            gadm_layer_id=config["build_shape_options"]["gadm_layer_id"],
            cluster_options=config["cluster_options"],
        input:
            network=GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}.nc"),
            country_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/country_shapes.geojson"
            ),
            regions_onshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_onshore_elec_s{simpl}.geojson"
            ),
            regions_offshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_offshore_elec_s{simpl}.geojson"
            ),
            #gadm_shapes="resources/" + RDIR + "shapes/MAR2.geojson",
            #using this line instead of the following will test updated gadm shapes for MA.
            #To use: downlaod file from the google drive and place it in resources/" + RDIR + "shapes/
            #Link: https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM
            gadm_shapes=GS.remote(
                BUCKET + "resources/" + RDIR + "shapes/gadm_shapes.geojson"
            ),
            # busmap=ancient('resources/" + RDIR + "bus_regions/busmap_elec_s{simpl}.csv'),
            # custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            #                if config["enable"].get("custom_busmap", False) else []),
            tech_costs=COSTS,
        output:
            network=GS.remote(
                BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}.nc"
            ),
            regions_onshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson"
            ),
            regions_offshore=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson"
            ),
            busmap=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/busmap_elec_s{simpl}_{clusters}.csv"
            ),
            linemap=GS.remote(
                BUCKET
                + "resources/"
                + RDIR
                + "bus_regions/linemap_elec_s{simpl}_{clusters}.csv"
            ),
        log:
            GS.remote(
                BUCKET
                + "logs/"
                + RDIR
                + "cluster_network/elec_s{simpl}_{clusters}.log"
            ),
        benchmark:
            GS.remote(
                BUCKET
                + "benchmarks/"
                + RDIR
                + "cluster_network/elec_s{simpl}_{clusters}"
            )
        threads: 1
        resources:
            mem_mb=3000,
        script:
            "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network=GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}.nc"),
        tech_costs=COSTS,
    output:
        GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec.nc"),
    log:
        GS.remote(
            BUCKET
            + "logs/"
            + RDIR
            + "add_extra_components/elec_s{simpl}_{clusters}.log"
        ),
    benchmark:
        GS.remote(
            BUCKET
            + "benchmarks/"
            + RDIR
            + "add_extra_components/elec_s{simpl}_{clusters}_ec"
        )
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_extra_components.py"


rule prepare_network:
    params:
        links=config["links"],
        lines=config["lines"],
        s_max_pu=config["lines"]["s_max_pu"],
        electricity=config["electricity"],
        costs=config["costs"],
    input:
        GS.remote(BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec.nc"),
        tech_costs=COSTS,
    output:
        GS.remote(
            BUCKET + "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
        ),
    log:
        GS.remote(
            BUCKET
            + "logs/"
            + RDIR
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.log"
        ),
    benchmark:
        GS.remote(
            BUCKET
            + "benchmarks/"
            + RDIR
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/prepare_network.py"


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters.endswith("flex"):
        return int(factor * (18000 + 180 * int(w.clusters[:-4])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    elif w.clusters == "min":
        return int(factor * (18000 + 180 * 20))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == False:

    rule solve_network:
        params:
            countries = config["countries"],
            solving=config["solving"],
            augmented_line_connection=config["augmented_line_connection"],
        input:
            network=GS.remote(
                BUCKET
                + "networks/"
                + RDIR
                + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
            ),
        output:
            network=GS.remote(
                BUCKET
                + "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc"
            ),
        log:
            solver=GS.remote(
                normpath(
                    BUCKET
                    + "logs/"
                    + RDIR
                    + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
                )
            ),
            python=GS.remote(
                BUCKET
                + "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log"
            ),
        benchmark:
            GS.remote(
                BUCKET
                + "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}"
            )
        threads: 20
        resources:
            mem=memory,
        shadow:
            "copy-minimal" if os.name == "nt" else "shallow"
        script:
            "scripts/solve_network.py"


if config["monte_carlo"]["options"].get("add_to_snakefile", False) == True:

    rule monte_carlo:
        params:
            monte_carlo=config["monte_carlo"],
        input:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        output:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        log:
            "logs/"
            + RDIR
            + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "prepare_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}"
            )
        threads: 1
        resources:
            mem_mb=4000,
        script:
            "scripts/monte_carlo.py"

    rule solve_monte:
        input:
            expand(
                "networks/"
                + RDIR
                + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
                **config["scenario"],
            ),

    rule solve_network:
        params:
            solving=config["solving"],
            augmented_line_connection=config["augmented_line_connection"],
        input:
            "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        output:
            "results/"
            + RDIR
            + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
        log:
            solver=normpath(
                "logs/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_solver.log"
            ),
            python="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_python.log",
            memory="logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}_memory.log",
        benchmark:
            (
                "benchmarks/"
                + RDIR
                + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}"
            )
        threads: 20
        resources:
            mem_mb=memory,
        shadow:
            "shallow"
        script:
            "scripts/solve_network.py"

    rule solve_all_networks_monte:
        input:
            expand(
                "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{unc}.nc",
                **config["scenario"],
            ),


def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return [COSTS] + expand(
        "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        ll=ll,
        **{
            k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
            for k in ["simpl", "clusters", "opts"]
        },
    )


rule make_summary:
    params:
        electricity=config["electricity"],
        costs=config["costs"],
        ll=config["scenario"]["ll"],
        scenario=config["scenario"],
    input:
        input_make_summary,
        tech_costs=COSTS,
    output:
        directory(
            "results/"
            + RDIR
            + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
        ),
    log:
        "logs/"
        + RDIR
        + "make_summary/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    script:
        "scripts/make_summary.py"


rule plot_summary:
    input:
        "results/"
        + RDIR
        + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
    output:
        "results/"
        + RDIR
        + "plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_summary/{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log",
    script:
        "scripts/plot_summary.py"


rule plot_network:
    params:
        electricity=config["electricity"],
        costs=config["costs"],
        plotting=config["plotting"],
    input:
        network="results/"
        + RDIR
        + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        africa_shape="resources/" + RDIR + "shapes/africa_shape.geojson",
        tech_costs=COSTS,
    output:
        only_map="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log",
    script:
        "scripts/plot_network.py"


rule build_test_configs:
    input:
        base_config="config.tutorial.yaml",
        update_file_list=[
            "test/config.tutorial_noprogress.yaml",
            "test/config.custom.yaml",
            "test/config.monte_carlo.yaml",
            "test/config.landlock.yaml",
        ],
    output:
        tmp_test_configs=[
            "test/tmp/config.tutorial_noprogress_tmp.yaml",
            "test/tmp/config.custom_tmp.yaml",
            "test/tmp/config.monte_carlo_tmp.yaml",
            "test/tmp/config.landlock_tmp.yaml",
        ],
    script:
        "scripts/build_test_configs.py"


rule make_statistics:
    params:
        countries=config["countries"],
        renewable_carriers=config["electricity"]["renewable_carriers"],
        renewable=config["renewable"],
        crs=config["crs"],
        scenario=config["scenario"],
    output:
        stats="results/" + RDIR + "stats.csv",
    threads: 1
    script:
        "scripts/make_statistics.py"


rule run_scenario:
    input:
        diff_config="configs/scenarios/config.{scenario_name}.yaml",
    output:
        touchfile=touch("results/{scenario_name}/scenario.done"),
        copyconfig="results/{scenario_name}/config.yaml",
    threads: 1
    resources:
        mem_mb=5000,
    run:
        from build_test_configs import create_test_config
        import yaml

        # get base configuration file from diff config
        with open(input.diff_config) as f:
            base_config_path = (
                yaml.full_load(f)
                .get("run", {})
                .get("base_config", "config.default.yaml")
            )

            # Ensure the scenario name matches the name of the configuration
        create_test_config(
            input.diff_config,
            {"run": {"name": wildcards.scenario_name}},
            input.diff_config,
        )
        # merge the default config file with the difference
        create_test_config(base_config_path, input.diff_config, "config.yaml")
        os.system("snakemake -j all solve_all_networks --rerun-incomplete")
        os.system("snakemake -j1 make_statistics --force")
        copyfile("config.yaml", output.copyconfig)



rule run_all_scenarios:
    input:
        expand(
            "results/{scenario_name}/scenario.done",
            scenario_name=[
                c.stem.replace("config.", "")
                for c in Path("configs/scenarios").glob("config.*.yaml")
            ],
        ),
rule plot_trace:
    params: 
        year = int(config["snapshots"]["start"].split("-")[0]),
        country = config["countries"]
    input:
        network=GS.remote(
                BUCKET
                + "results/"
                + RDIR
                + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc"
            ),
        ember_annual = GS.remote("feo-pypsa-staging/data/GER2024-Data-Yearly.csv")
    output:
        annual_generation=GS.remote(
                BUCKET
                + "results/"
                + RDIR
                + "plots/annual_generation_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.png"

            ),
    log:
        GS.remote(
            BUCKET
            + "logs/"
            + RDIR
            + "plot_trace/annual_generation{simpl}_{clusters}_ec_l{ll}_{opts}_trace.log"
        ),
    script:
        "scripts/plot_trace.py"

rule plot_all_trace:
    input:
        GS.remote(expand(
            BUCKET
            + "results/"
            + RDIR
            + "plots/annual_generation_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.png",
            **config["scenario"],
        )),