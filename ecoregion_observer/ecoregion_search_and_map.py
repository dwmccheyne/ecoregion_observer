#!/usr/bin/env python3

import os
import requests
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import zipfile
import argparse
import time
import urllib3
import matplotlib.pyplot as plt

# Suppress SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Constants
INATURALIST_API_URL = "https://api.inaturalist.org/v1/observations"
EPA_BASE_URL = "https://gaftp.epa.gov/epadatacommons/ORD/Ecoregions"
DEFAULT_SHAPEFILE_DIR = os.getenv("SHAPEFILE_DIR", "./data/ecoregions")


def fetch_shapefiles(states, destination, level):
    """
    Download and extract ecoregion shapefiles for specified states from the EPA website.
    """
    os.makedirs(destination, exist_ok=True)
    for state in states:
        state_lower = state.lower()
        level_suffix = f"_eco_{level.lower()}"
        url = f"{EPA_BASE_URL}/{state_lower}/{state_lower}{level_suffix}.zip"
        zip_path = os.path.join(destination, f"{state_lower}{level_suffix}.zip")
        try:
            response = requests.get(url, stream=True, verify=False)
            if response.status_code == 200:
                with open(zip_path, "wb") as f:
                    f.write(response.content)
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(destination)
                print(f"Downloaded and extracted data for {state} (Level {level}).")
            else:
                print(f"Failed to fetch data for {state} (Status Code: {response.status_code}).")
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data for {state}: {e}")


def load_shapefiles_from_states(states, directory, level):
    """
    Load and combine shapefiles for specified states and ecoregion level.
    """
    level_suffix = f"_eco_{level.lower()}.shp"
    shapefiles = [
        os.path.join(directory, f"{state.lower()}{level_suffix}") for state in states
    ]
    for shp in shapefiles:
        if not os.path.exists(shp):
            raise FileNotFoundError(f"Shapefile not found: {shp}")
    gdfs = []
    for shp in shapefiles:
        gdf = gpd.read_file(shp)
        if gdf.crs != "EPSG:4326":
            gdf = gdf.to_crs("EPSG:4326")
        gdfs.append(gdf)
    return gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True))


def get_bounding_box_from_shapefile(shapefile):
    """
    Extract the bounding box coordinates from a given shapefile.
    """
    gdf = gpd.read_file(shapefile)
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = gdf.total_bounds
    return {"swlat": miny, "swlng": minx, "nelat": maxy, "nelng": maxx}


def fetch_observations(taxon_id, filters, per_page=200, max_results=10000, bounding_box=None, place_id=None):
    """
    Fetch observations from the iNaturalist API based on specified criteria.
    """
    observations = []
    page = 1
    while len(observations) < max_results:
        params = {
            "taxon_id": taxon_id,
            "per_page": per_page,
            "page": page,
            "geo": "true",
        }
        for f_key, f_value in filters.items():
            params[f_key] = f_value
        if place_id:
            params["place_id"] = place_id
        elif bounding_box:
            params.update(bounding_box)
        response = requests.get(INATURALIST_API_URL, params=params)
        if response.status_code != 200:
            raise requests.exceptions.RequestException(
                f"Failed to fetch data from iNaturalist API (status code {response.status_code})."
            )
        data = response.json()
        results = data.get("results", [])
        observations.extend(results)
        print(f"Page {page}: Fetched {len(results)} observations. Total so far: {len(observations)}")
        if len(results) < per_page:
            break
        page += 1
        time.sleep(0.1)
    return observations[:max_results]


def search_observations(taxon_id, states, filters, shp_dir, output_file, level, place_id=None):
    """
    Searches for observations of a given taxon within specified states and ecoregions.
    Saves results to a CSV file.
    """
    ecoregions = load_shapefiles_from_states(states, shp_dir, level)
    shapefile = os.path.join(shp_dir, f"{states[0].lower()}_eco_{level.lower()}.shp")
    bounding_box = get_bounding_box_from_shapefile(shapefile)
    observations = fetch_observations(
        taxon_id, filters, bounding_box=bounding_box, place_id=place_id
    )
    obs_data = {
        "id": [obs["id"] for obs in observations],
        "geometry": [Point(obs["geojson"]["coordinates"]) for obs in observations if "geojson" in obs],
        "scientific_name": [obs.get("taxon", {}).get("name", None) for obs in observations],
    }
    obs_gdf = gpd.GeoDataFrame(obs_data, crs="EPSG:4326")
    if ecoregions.crs != obs_gdf.crs:
        ecoregions = ecoregions.to_crs(obs_gdf.crs)
    joined = gpd.sjoin(obs_gdf, ecoregions, how="left", predicate="within")
    ecoregion_name_column = next((col for col in ecoregions.columns if "NAME" in col.upper()), None)
    if not ecoregion_name_column:
        raise KeyError("No ecoregion name column found in the ecoregions data.")
    group = joined.groupby(ecoregion_name_column)
    counts = group.size().reset_index(name="observation_count")
    unique_species = group["scientific_name"].agg(lambda x: x.nunique()).reset_index(name="unique_species_count")
    all_ecoregions = ecoregions[[ecoregion_name_column]].drop_duplicates()
    full_counts = all_ecoregions.merge(counts, on=ecoregion_name_column, how="left")
    full_counts = full_counts.merge(unique_species, on=ecoregion_name_column, how="left")
    full_counts["observation_count"] = full_counts["observation_count"].fillna(0).astype(int)
    full_counts["unique_species_count"] = full_counts["unique_species_count"].fillna(0).astype(int)
    full_counts.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}.")


def visualize_ecoregion_data(csv_file, shapefile_dir, output_png, level):
    """
    Visualize observation data on a map using ecoregion shapefiles.
    """
    try:
        data = pd.read_csv(csv_file)
        print(f"Loaded data: {data.head()}")
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {csv_file}")

    level_suffix = f"_eco_{level.lower()}.shp"
    shapefiles = [
        os.path.join(shapefile_dir, fname)
        for fname in os.listdir(shapefile_dir)
        if fname.endswith(level_suffix)
    ]
    if not shapefiles:
        raise FileNotFoundError(f"No shapefiles found for level {level} in {shapefile_dir}.")
    
    ecoregions = gpd.GeoDataFrame(
        pd.concat([gpd.read_file(shp) for shp in shapefiles], ignore_index=True)
    )
    if ecoregions.crs is None:
        ecoregions = ecoregions.set_crs("EPSG:4326")
    elif ecoregions.crs != "EPSG:4326":
        ecoregions = ecoregions.to_crs("EPSG:4326")

    ecoregion_name_column = next((col for col in ecoregions.columns if "NAME" in col.upper()), None)
    if not ecoregion_name_column:
        raise KeyError("No ecoregion name column found in the shapefiles.")

    merged = ecoregions.merge(data, left_on=ecoregion_name_column, right_on=ecoregion_name_column, how="left")
    merged["observation_count"] = merged["observation_count"].fillna(0).astype(int)
    merged["unique_species_count"] = merged["unique_species_count"].fillna(0).astype(int)

    fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    merged.boundary.plot(ax=ax[0], linewidth=0.5, color="black")
    merged.plot(column="observation_count", cmap="Oranges", legend=True, ax=ax[0], edgecolor="black")
    ax[0].set_title("Observation Counts by Ecoregion")

    merged.boundary.plot(ax=ax[1], linewidth=0.5, color="black")
    merged.plot(column="unique_species_count", cmap="Purples", legend=True, ax=ax[1], edgecolor="black")
    ax[1].set_title("Unique Species Counts by Ecoregion")

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"Visualization saved as {output_png}")


def display_help():
    """
    Print a detailed help message for the script.
    """
    print("""
    Usage: ./script_name.py <command> [options]

    Commands:
        fetch       Download and extract ecoregion shapefiles.
        search      Search for iNaturalist observations in specified ecoregions.
        visualize   Create map visualizations of observation data.

    Use the -h or --help option with any command for details.
    """)


def main():
    parser = argparse.ArgumentParser(description="Analyze iNaturalist observations by ecoregion.")
    subparsers = parser.add_subparsers(dest="command")

    fetch_parser = subparsers.add_parser("fetch", help="Fetch ecoregion shapefiles for specified states.")
    fetch_parser.add_argument("states", nargs="+", help="List of state abbreviations.")
    fetch_parser.add_argument("--output", default=DEFAULT_SHAPEFILE_DIR, help="Destination directory for shapefiles.")
    fetch_parser.add_argument("--level", choices=["L3", "L4"], default="L3", help="Ecoregion level (L3 or L4).")

    search_parser = subparsers.add_parser("search", help="Search iNaturalist observations by ecoregion.")
    search_parser.add_argument("taxon_id", type=int, help="Taxon ID to query.")
    search_parser.add_argument("--states", nargs="+", required=True, help="List of states to include in the analysis.")
    search_parser.add_argument("--filter", action="append", help="Additional filters for the API query.")
    search_parser.add_argument("--shp-dir", default=DEFAULT_SHAPEFILE_DIR, help="Directory containing shapefiles.")
    search_parser.add_argument("--level", choices=["L3", "L4"], default="L3", help="Ecoregion level (L3 or L4).")
    search_parser.add_argument("--output", default="observations_by_ecoregion.csv", help="Output file.")

    visualize_parser = subparsers.add_parser("visualize", help="Visualize observation data on a map.")
    visualize_parser.add_argument("csv_file", help="Path to CSV file with observation data.")
    visualize_parser.add_argument("shapefile_dir", help="Directory containing ecoregion shapefiles.")
    visualize_parser.add_argument("output_png", help="Output PNG file for the visualization.")
    visualize_parser.add_argument("--level", choices=["L3", "L4"], default="L3", help="Ecoregion level (L3 or L4).")

    args = parser.parse_args()

    if args.command == "fetch":
        fetch_shapefiles(args.states, args.output, args.level)
    elif args.command == "search":
        filters = {}
        if args.filter:
            for f in args.filter:
                key, value = f.split("=", 1)
                filters[key.strip()] = value.strip()
        search_observations(
            taxon_id=args.taxon_id,
            states=args.states,
            filters=filters,
            shp_dir=args.shp_dir,
            output_file=args.output,
            level=args.level,
        )
    elif args.command == "visualize":
        visualize_ecoregion_data(args.csv_file, args.shapefile_dir, args.output_png, args.level)
    else:
        display_help()


if __name__ == "__main__":
    main()
