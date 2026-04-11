import sys
import os

import matplotlib
from eoio import read

import s2_rut_python._vendor  # noqa: F401

from S2RUT import S2RUT_L1

import os

this_directory = os.path.dirname(__file__)

file = os.path.abspath(
    os.path.join(
        os.path.dirname(this_directory),
        "third-party",
        "Data",
        "S2A_MSIL1C_20210310T084801_N0500_R107_T33KWP_20230525T152740.SAFE",
    )
)
# Open product using eoio
ds = read(
    file,
    vars_sel={
        "meas": ["B01"],  # Read only a few bands for efficiency
        "aux": [
            "solar_zenith_angle",
            "solar_azimuth_angle",
            "viewing_zenith_angle",
            "viewing_azimuth_angle",
        ],
    },
    read_params={
        "use_chunks": False,
        "metadata_level": "all",
        "save_extracted": True,
        "ave_va_det": True,  # Average over the detector dimension for the angle variables
    },
    processors={
        "interpolate": {
            "coords": ["y_5000m", "x_5000m"],
            "data_vars": ["solar_zenith_angle"],
            "target_grid": ["y_60m", "x_60m"],
            "method": "linear",
            "inplace": True,
        },
    },
)

print(ds.data_vars)
# print(ds['solar_zenith_angle'])
print(ds.coords)

# Compute uncertainty (absolute value, in reflectance dimension)
input_contributors = os.path.abspath(
    os.path.join(
        os.path.dirname(this_directory), "third-party", "Data", "unc_contributors.json"
    )
)

metadata = {
    "spacecraft": None,
    "quant": None,
    "A": [],
    "offset": [],
    "alpha": {},
    "beta": {},
    "Esun": [],
    "Usun": None,
    "refined": False,
}  # this dictionary contains the relevant metadata parameters
# read required metadata from ds.attrs and fill the metadata dictionary
METADATA_MAP = {
    "spacecraft": "platform",
    "quant": "quantification_level",
    "Usun": "reflectance_conversion_u",
}
VAR_MTD_MAP = {
    "Esun": "solar_irradiance",
    "alpha": "noise_model_alpha",
    "beta": "noise_model_beta",
    "A": "physical_gains",
    "offset": "radiometric_offset",
}

for key, path in METADATA_MAP.items():
    if path in ds.attrs:
        metadata[key] = ds.attrs[path]
    elif path in ds.attrs["product_metadata"]:
        metadata[key] = ds.attrs["product_metadata"][path]
    else:
        raise KeyError(
            f"Required metadata parameter '{key}' not found in dataset attributes at path '{path}'."
        )

for key, path in VAR_MTD_MAP.items():
    metadata[key] = {
        band: ds[band].attrs["product_metadata"][path]
        for band in ds.data_vars
        if "B" in band
    }

print("Metadata extracted for uncertainty calculation:", metadata)
RUTl1 = S2RUT_L1(input_contributors)
u_ref, u_cont = RUTl1.unc_calculation_abs(
    ds["B01"].values,
    "B01",
    0,
    metadata,
    [ds["solar_zenith_angle"].values],
    do_contributor=True,
)


print("Uncertainty calculated for B01:", u_ref, u_cont)
print("Shape of u_ref:", u_ref.shape)
print("u_cont keys:", u_cont.keys())


from s2_rut_python.interface import S2RUTTool

rut_processor = S2RUTTool()
ds = rut_processor.run(
    ds,
    data_vars=["B01"],
    group_unc=False,  # set to True to compute group uncertainties (systematic and random)
)

# Check the resulting dataset
print("Data variables in resulting dataset:", ds.data_vars)
print(ds.unc)

# compare with reference uncertainty
import matplotlib.pyplot as plt

# plt.imshow(ds.B01.values)
# plt.title("B01 reflectance")
# plt.colorbar()
# plt.show()
# for unc in ds.data_vars:
#     if 'u_' in unc:
#         print(f"Comparing {unc} with reference uncertainty...")
#         if unc.split('_B01')[0] in u_cont:
#             plt.imshow(ds[unc].values - u_cont[unc.split('_B01')[0]])
#             plt.title(f"Difference between {unc} and reference uncertainty")
#             plt.colorbar()
#             plt.show()

ds = read(
    file,
    vars_sel={
        "meas": ["B01"],  # Read only a few bands for efficiency
        "aux": [
            "solar_zenith_angle",
            "solar_azimuth_angle",
            "viewing_zenith_angle",
            "viewing_azimuth_angle",
        ],
    },
    read_params={
        "use_chunks": False,
        "metadata_level": "all",
        "save_extracted": True,
        "ave_va_det": True,  # Average over the detector dimension for the angle variables
    },
    processors={
        "interpolate": {
            "coords": ["y_5000m", "x_5000m"],
            "data_vars": ["solar_zenith_angle"],
            "target_grid": ["y_60m", "x_60m"],
            "method": "linear",
            "inplace": True,
        },
        "s2_rut": {
            "data_vars": ["B01"],
            "group_unc": True,  # set to True to compute group uncertainties (systematic and random)
        },
    },
)

print("Data variables in resulting dataset after s2_rut processor:", ds.data_vars)
print(ds.unc)

for unc in ds.data_vars:
    if "u_" in unc:
        print(f"Comparing {unc} with reference uncertainty...")
        if unc.split("_B01")[0] in u_cont:
            plt.imshow(ds[unc].values - u_cont[unc.split("_B01")[0]])
            plt.title(f"Difference between {unc} and reference uncertainty")
            plt.colorbar()
            plt.show()
        else:
            import numpy as np

            (100 * ds[unc] / ds["B01"]).plot(
                robust=True
            )  # , vmin=np.nanpercentile(100*ds[unc].values / ds['B01'].values, 1), vmax=np.nanpercentile(100*ds[unc].values / ds['B01'].values, 99))  # relative uncertainty in percentage
            plt.title(f"{unc} (no reference uncertainty available)")
            # plt.colorbar()
            plt.show()
