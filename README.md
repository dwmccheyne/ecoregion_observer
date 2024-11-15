# Ecoregion Observer

**Ecoregion Observer** is a Python package for fetching, analyzing, and visualizing ecoregion-based observations using data from the iNaturalist API and EPA's ecoregion shapefiles.

## Features

- Fetch EPA ecoregion shapefiles for specified states and levels.
- Query iNaturalist API for observations within those ecoregions.
- Visualize observation counts and unique species counts on a map.

## Installation

1. Clone the repository or install the package via pip:
   ```bash
   pip install git+https://github.com/dwmccheyne/ecoregion_observer.git
   ```

2. Ensure required Python packages are installed (automatically handled by pip):
   - `requests`
   - `geopandas`
   - `pandas`
   - `matplotlib`
   - `shapely`

## Command-Line Usage

The package provides three primary commands: `fetch`, `search`, and `visualize`.

### Example Workflow

Below is a complete example using the package to fetch ecoregion shapefiles for New York, search for fungi observations, and visualize the results.

### Step 1: Fetch Ecoregion Shapefiles
Download Level IV (L4) ecoregion shapefiles for New York:
```bash
ecoregion-search-and-map fetch NY --output ./data/ecoregions --level L4
```

### Step 2: Search for Observations
Search for fungi observations in New York's ecoregions with a DNA Barcode ITS filter and save the results to `fungi_observations.csv`:
```bash
ecoregion-search-and-map search 47170 --states NY --shp-dir ./data/ecoregions --output fungi_observations.csv --level L4 --filter "field:DNA Barcode ITS="
```
- `47170` is the taxon ID for fungi. Replace it with other IDs for different taxa if needed.
- The `--filter` option can include additional iNaturalist filters.

### Step 3: Visualize the Results
Create a visualization of fungi observations, showing observation counts and unique species counts across New York's ecoregions:
```bash
ecoregion-search-and-map visualize fungi_observations.csv ./data/ecoregions fungi_visualization.png --level L4
```
The output map will be saved as `fungi_visualization.png`.

## Example Python Usage

In addition to command-line use, the package can be used programmatically:

```python
from ecoregion_search_and_map.ecoregion_search_and_map import fetch_shapefiles, search_observations, visualize_ecoregion_data

# Step 1: Fetch shapefiles
fetch_shapefiles(states=["NY"], destination="./data/ecoregions", level="L4")

# Step 2: Search for observations
filters = {"field:DNA Barcode ITS": "yes"}
search_observations(
    taxon_id=47170,
    states=["NY"],
    filters=filters,
    shp_dir="./data/ecoregions",
    output_file="fungi_observations.csv",
    level="L4"
)

# Step 3: Visualize the results
visualize_ecoregion_data(
    csv_file="fungi_observations.csv",
    shapefile_dir="./data/ecoregions",
    output_png="fungi_visualization.png",
    level="L4"
)
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.