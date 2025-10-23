""" Utilities functions for processors"""

import numpy as np
from scipy.interpolate import griddata  # type: ignore[import-untyped]


def interp_sza_s2(ds, target_resolution, intrp_method="linear"):
    """
    Interpolate the SZA from 5000m to target resolution

    :param ds: S2 ds read in by eoio
    :param target_resolution: resolution to interpolate to
    :return: interpolated SZA on target resolution grid
    """
    # check angles are in the dataset
    if "x_5000m" not in ds.dims or "y_5000m" not in ds.dims:
        raise KeyError(
            "S2 Radiance Processor requires observation geometry to be read in - check this is included in aux_data."
        )
    # get source x and y
    x_source = ds["x_5000m"].data
    y_source = ds["y_5000m"].data

    # get target x and y
    if (
        f"x_{target_resolution:.0f}m" not in ds.dims
        or f"y_{target_resolution:.0f}m" not in ds.dims
    ):
        raise KeyError(
            f"Resolution {target_resolution} not included in dataset. Existing dimensions are {ds.dims}"
        )
    x_output = ds[f"x_{target_resolution:.0f}m"].data
    y_output = ds[f"y_{target_resolution:.0f}m"].data

    # create source and target grids in format needed for griddata
    x_source_grid, y_source_grid = np.meshgrid(x_source, y_source)
    x_output_grid, y_output_grid = np.meshgrid(x_output, y_output)

    # get solar zenith angle data from dataset
    if "solar_zenith_angle" in ds:
        data = ds["solar_zenith_angle"].data
    elif "Zenith_Sun_Angles_Grid" in ds:
        data = ds["Zenith_Sun_Angles_Grid"].data
    else:
        raise KeyError(
            f"No variable named 'Zenith_Sun_Angles' or 'Zenith_Sun_Angles_Grid'. Variables on the dataset include {ds.data_vars}"
        )

    # interpolate values using source and desired x y output
    output_data = griddata(
        (x_source_grid.flatten(), y_source_grid.flatten()),
        data.flatten(),
        (x_output_grid, y_output_grid),
        method=intrp_method,
    )

    return output_data


if __name__ == "__main__":
    pass
