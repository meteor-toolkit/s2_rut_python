"""Utility helpers for examples."""

import os
import re

import numpy as np
from eoio import read


def infer_roi_crs(ds_probe, probe_band: str):
    """Try to infer ROI CRS from band/grid-mapping metadata."""
    roi_crs = ds_probe[probe_band].attrs.get("crs") or ds_probe.attrs.get("crs")
    if roi_crs:
        return roi_crs

    grid_mapping = ds_probe[probe_band].attrs.get("grid_mapping")
    if grid_mapping in ds_probe:
        gm_attrs = ds_probe[grid_mapping].attrs

        if "epsg_code" in gm_attrs:
            return f"EPSG:{gm_attrs['epsg_code']}"

        if "spatial_ref" in gm_attrs:
            wkt = gm_attrs["spatial_ref"]
            match = re.search(r"EPSG\"\s*,\s*\"(\d+)\"", wkt)
            if match:
                return f"EPSG:{match.group(1)}"

    return None


def infer_roi_crs_from_product_name(safe_product: str):
    """Infer UTM EPSG from Sentinel-2 tile ID in product name (e.g. T31TFJ)."""
    product_name = os.path.basename(safe_product)
    match = re.search(r"_T(\d{2})([A-Z])[A-Z]{2}_", product_name)
    if not match:
        return None

    zone = int(match.group(1))
    lat_band = match.group(2)
    hemisphere_north = lat_band >= "N"
    epsg = 32600 + zone if hemisphere_north else 32700 + zone
    return f"EPSG:{epsg}"


def build_centered_roi_subset(
    safe_product: str, probe_band: str = "B01", roi_size_m: float = 1000.0
):
    """Build a centered square ROI subset dictionary (~roi_size_m by roi_size_m)."""
    ds_probe = read(
        safe_product,
        vars_sel={"meas": [probe_band]},
        read_params={
            "use_chunks": True,
            "metadata_level": "all",
            "save_extracted": True,
            "ave_va_det": True,
        },
    )

    y_dim = next(dim for dim in ds_probe[probe_band].dims if dim.startswith("y_"))
    x_dim = next(dim for dim in ds_probe[probe_band].dims if dim.startswith("x_"))

    x_vals = ds_probe.coords[x_dim].values
    y_vals = ds_probe.coords[y_dim].values
    x_min, x_max = float(np.min(x_vals)), float(np.max(x_vals))
    y_min, y_max = float(np.min(y_vals)), float(np.max(y_vals))

    x_center = 0.5 * (x_min + x_max)
    y_center = 0.5 * (y_min + y_max)
    half = 0.5 * roi_size_m

    roi = (
        max(x_min, x_center - half),
        max(y_min, y_center - half),
        min(x_max, x_center + half),
        min(y_max, y_center + half),
    )

    subset = {"roi": roi}
    roi_crs = infer_roi_crs(ds_probe, probe_band)
    if roi_crs is None:
        roi_crs = infer_roi_crs_from_product_name(safe_product)
    if roi_crs:
        subset["roi_crs"] = roi_crs

    return subset
