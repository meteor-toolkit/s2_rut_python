"""example: read a Sentinel-2 L1C SAFE product and run S2-RUT."""

import os
import glob
import socket

import matplotlib.pyplot as plt
from eoio import read
from shapely import wkt

from s2_rut_python.interface import S2RUTTool

if (
    socket.gethostname() == "lyon.npl.co.uk"
    or socket.gethostname() == "leipzig.npl.co.uk"
    or "leiden" in socket.gethostname()
):
    DATA_DIRECTORY = "/mnt/t/data/"
else:
    DATA_DIRECTORY = r"T:\ECO\EOServer\data"


# Define input path from test datasets
SAFE_PATH = os.path.join(
    DATA_DIRECTORY,
    "unittest_datasets",
    "S2MSIL1C",
    "S2A_MSIL1C_20251128T111431_N0511_R137_T30UXE_20251128T121631.SAFE",
)
print("SAFE PATH:", SAFE_PATH)

# Define ROI
roi_string = "POLYGON ((\
 -0.913910 53.786950,\
 -0.382900 53.776490,\
 -0.406240 53.464690,\
 -0.925060 53.474770,\
 -0.913910 53.786950\
))"

geom = wkt.loads(roi_string)

# Read the dataset with eoio, selecting only a few bands and auxiliary variables for efficiency.
bands = ["B01", "B03", "B09"]
ds = read(
    SAFE_PATH,
    vars_sel={
        "meas": bands,
        "aux": [
            "solar_zenith_angle",
            "solar_azimuth_angle",
            "viewing_zenith_angle",
            "viewing_azimuth_angle",
        ],
    },
    read_params={
        "use_chunks": True,
        "metadata_level": "all",
        "save_extracted": True,
        "ave_va_det": True,
    },
    subset={"roi": geom, "roi_crs": "EPSG:4326"},
    processors={
        "interpolate": {
            "coords": ["y_5000m", "x_5000m"],
            "data_vars": ["solar_zenith_angle"],
            "target_grid": [["y_60m", "y_10m"], ["x_60m", "x_10m"]],
            "method": "linear",
            "inplace": False,
        }
    },
)

# Run S2-RUT on the dataset, computing both random and systematic uncertainties.
print("Running S2-RUT")
rut = S2RUTTool()
ds_out = rut.run(ds=ds, data_vars=bands, group_unc=True)


# Create a single grid plot: rows = bands, columns = [reflectance, random [%], systematic [%]].
fig, axes = plt.subplots(
    nrows=len(bands),
    ncols=3,
    figsize=(15, 4 * len(bands)),
)

# Handle single band case (axes is 1D instead of 2D).
if len(bands) == 1:
    axes = axes.reshape(1, -1)

for i, band in enumerate(bands):
    # Column 0: Reflectance
    ds_out[band].plot(ax=axes[i, 0], robust=True)
    axes[i, 0].set_title(f"{band} Reflectance")
    axes[i, 0].set_xlabel("")
    axes[i, 0].set_ylabel("")

    # Column 1: Random uncertainty [%]
    random_unc_pct = 100.0 * ds_out[f"u_random_{band}"] / ds_out[band]
    random_unc_pct.plot(ax=axes[i, 1], robust=True)
    axes[i, 1].set_title(f"{band} Random Uncertainty [%]")
    axes[i, 1].set_xlabel("")
    axes[i, 1].set_ylabel("")

    # Column 2: Systematic uncertainty [%]
    sys_unc_pct = 100.0 * ds_out[f"u_systematic_{band}"] / ds_out[band]
    sys_unc_pct.plot(ax=axes[i, 2], robust=True)
    axes[i, 2].set_title(f"{band} Systematic Uncertainty [%]")
    axes[i, 2].set_xlabel("")
    axes[i, 2].set_ylabel("")

plt.tight_layout()
plt.savefig("s2rut_example.png")


if __name__ == "__main__":
    pass
