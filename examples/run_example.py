"""Simple example: read a Sentinel-2 L1C SAFE product and run S2-RUT."""

import os
import glob
import socket

import matplotlib.pyplot as plt
from eoio import read

from example_utils import build_centered_roi_subset
from s2_rut_python.interface import S2RUTTool

if socket.gethostname() == "lyon.npl.co.uk" or socket.gethostname() == "leipzig.npl.co.uk" or 'leiden' in socket.gethostname():
    S2_DIRECTORY = "/mnt/t/data/product_archive/S2A/S2_MSI_L1C"
else:
    S2_DIRECTORY = r"T:\ECO\EOServer\data\product_archive\S2A\S2_MSI_L1C"

# Use T drive example products.
safe_candidates = glob.glob(
    os.path.join(
        S2_DIRECTORY,
        '*',
        '*',
        '*',
        "*.SAFE",
    ))

if len(safe_candidates) == 0:
    raise FileNotFoundError(
        "Could not find the example SAFE product on the T drive. "
        "Update `safe_candidates` to a valid local path."
    )


# Process multiple bands for richer example.
bands = ["B01", "B03", "B09"]

safe_product = safe_candidates[0]
print(f"Reading SAFE product: {safe_product}")

subset = build_centered_roi_subset(safe_product, probe_band=bands[0], roi_size_m=2000.0)
print(f"Using automated ROI subset: {subset}")


ds = read(
    safe_product,
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
    subset=subset,
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

print("Running S2-RUT")
rut = S2RUTTool()
ds_out = rut.run(ds=ds, data_vars=bands, group_unc=True)

print("Output variables:")
print(list(ds_out.data_vars))

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
