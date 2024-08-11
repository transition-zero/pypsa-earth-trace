import argparse

import yaml
from google.cloud import storage
from snakemake.io import expand

parser = argparse.ArgumentParser()
parser.add_argument("run", nargs="+", type=str)
parser.add_argument("--scenario", type=lambda x: yaml.safe_load(x), required=True)
args = parser.parse_args()

BUCKET = "feo-pypsa-staging"
ISO2 = {
    "landlocked": [
        "AF", "AM", "BF", "BI", "BO", "BT", "BW", "BY", "CF", "CM", "ET", 
        "GA", "IQ", "KG", "LA", "LS", "MD", "ML", "MN", "MW", "NE", "NP", 
        "PY", "RW", "SS", "SZ", "TD", "TJ", "UG", "UZ", "XK", "ZM", "ZW",
    ],
    "coastal": [
        "AE", "AG", "AO", "AR", "AS", "AU", "AW", "AZ", "BB", "BD", "BH", 
        "BJ", "BM", "BN", "BR", "BS", "BZ", "CA", "CD", "CG", "CI", "CK", 
        "CL", "CN", "CO", "CR", "CU", "CV", "CY", "DJ", "DM", "DO", "DZ", 
        "EC", "EG", "ER", "FJ", "FK", "FO", "GD", "GE", "GF", "GH", "GI", 
        "GL", "GM", "GN", "GP", "GQ", "GT", "GU", "GW", "GY", "HN", "HT", 
        "ID", "IL", "IN", "IR", "IS", "JM", "JO", "JP", "KE", "KH", "KI", 
        "KM", "KN", "KP", "KR", "KW", "KY", "KZ", "LB", "LC", "LK", "LR", 
        "LY", "MA", "MG", "MM", "MQ", "MR", "MS", "MT", "MU", "MV", "MX", 
        "MY", "MZ", "NA", "NC", "NG", "NI", "NR", "NU", "NZ", "OM", "PA", 
        "PE", "PF", "PG", "PH", "PK", "PR", "PS", "QA", "RE", "RU", "SA", 
        "SB", "SC", "SD", "SG", "SL", "SN", "SO", "SR", "ST", "SV", "SY", 
        "TC", "TG", "TH", "TL", "TM", "TN", "TO", "TR", "TT", "TW", "TZ", 
        "UA", "US", "UY", "VC", "VE", "VG", "VI", "VN", "VU", "WS", "YE", 
        "ZA",
    ],
}  # fmt: skip
FILES = [
    "cutouts/{run}/{cutout}.nc",
    "networks/{run}/base.nc",
    "networks/{run}/elec.nc",
    "networks/{run}/elec_s{simpl}.nc",
    "networks/{run}/elec_s{simpl}_{clusters}.nc",
    "networks/{run}/elec_s{simpl}_{clusters}_ec.nc",
    "networks/{run}/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    "resources/{run}/costs.csv",
    "resources/{run}/costs_modified.csv",
    "resources/{run}/demand_profiles.csv",
    "resources/{run}/natura.tiff",
    "resources/{run}/powerplants.csv",
    "resources/{run}/powerplants_osm2pm.csv",
    "resources/{run}/base_network/all_buses_build_network.csv",
    "resources/{run}/base_network/all_converters_build_network.csv",
    "resources/{run}/base_network/all_lines_build_network.csv",
    "resources/{run}/base_network/all_transformers_build_network.csv",
    "resources/{run}/bus_regions/busmap_elec_s{simpl}.csv",
    "resources/{run}/bus_regions/busmap_elec_s{simpl}_{clusters}.csv",
    "resources/{run}/bus_regions/connection_costs_s{simpl}.csv",
    "resources/{run}/bus_regions/linemap_elec_s{simpl}_{clusters}.csv",
    "resources/{run}/bus_regions/regions_offshore.geojson",
    "resources/{run}/bus_regions/regions_offshore_elec_s{simpl}.geojson",
    "resources/{run}/bus_regions/regions_offshore_elec_s{simpl}_{clusters}.geojson",
    "resources/{run}/bus_regions/regions_onshore.geojson",
    "resources/{run}/bus_regions/regions_onshore_elec_s{simpl}.geojson",
    "resources/{run}/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson",
    "resources/{run}/osm/clean/all_clean_generators.csv",
    "resources/{run}/osm/clean/all_clean_generators.geojson",
    "resources/{run}/osm/clean/all_clean_lines.geojson",
    "resources/{run}/osm/clean/all_clean_substations.geojson",
    "resources/{run}/osm/raw/all_raw_cables.geojson",
    "resources/{run}/osm/raw/all_raw_generators.csv",
    "resources/{run}/osm/raw/all_raw_generators.geojson",
    "resources/{run}/osm/raw/all_raw_lines.geojson",
    "resources/{run}/osm/raw/all_raw_substations.geojson",
    "resources/{run}/renewable_profiles/profile_{technology}.nc",
    "resources/{run}/shapes/africa_shape.geojson",
    "resources/{run}/shapes/country_shapes.geojson",
    "resources/{run}/shapes/gadm_shapes.geojson",
    "resources/{run}/shapes/offshore_shapes.geojson",
    "results/{run}/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_trace.nc",
]

storage_client = storage.Client()
bucket = storage_client.get_bucket(BUCKET)
blobs = [
    blob.name
    for run in args.run
    # NOTE: [!l!b] excludes the "logs" and "benchmarks" prefixes
    for blob in bucket.list_blobs(match_glob=f"[!l!b]*/{run}/**")
]
expected_blobs = expand(
    FILES,
    run=args.run,
    cutout=["cutout-era5"],
    technology=["hydro", "offwind-ac", "offwind-dc", "onwind", "solar"],
    **args.scenario,
)
missing_blobs = set(expected_blobs) - set(blobs)
other_blobs = set(blobs) - set(expected_blobs)

print("Missing blobs:", " ".join(missing_blobs))
print("Other blobs:", " ".join(other_blobs))
