"""
Code utilising eoio to read S2 L1 data product and apply s2-rut algorithm to calculate uncertainties of the product.

Data set selected is from the eoio documentation, read_examples_s2.py

Created 02/05/2024 by R.Ormane
"""
import os
import socket
import datetime

import matplotlib.pyplot as plt

import s2_rut_algo
import s2_l1_rad_conf
from eoio.interface import (
    read,
    product_bounds,
    product_processors,
    product_subsetting_params,
    mid_lon_lat,
)
from eoio.utils.dict_tools import *
from eoio.processors import utils as util

__all__ = []

s2_l1c_filepaths = os.path.abspath(r"T:\ECO\EOServer\data\satellite\S2A_MSI\L1\PICS\Libya4\S2A_MSIL1C_20190113T090331_N0207_R007_T34RGS_20190113T111955.SAFE")

s2_l1c_product_bounds = product_bounds(s2_l1c_filepaths)
s2_l1c_mid_lat_lon = mid_lon_lat(s2_l1c_filepaths)

s2_l1c_ds = read(
    s2_l1c_filepaths,
    subset_info={
        "meas_vars": ["B01"],
        "metadataLevel": True,
        "masks": ["ancillary_lost", "opaque", "cirrus", "detector_footprint"],
        "aux_data": [
            "observation_geometry",
            "tcwv",
            "msl",
        ],
    },
    process_params={
        "convert": ["reflectance"],
        "angles": {"angle_type": 'both'},
    },
)


product = 'Sentinel-2A'
band = 'B01'
subset = s2_l1c_ds['B01'].values  # [1000:1501, 1000:1501]

# Data from Variable Metadata

a = s2_l1c_ds[band].PHYSICAL_GAINS  # physical gains
e_sun = s2_l1c_ds[band].SOLAR_IRRADIANCE['#text']  # solar irradiance
alpha = s2_l1c_ds[band].ALPHA  # noise
beta = s2_l1c_ds[band].BETA  # noise

# Data from Global Metadata

u_sun = get_value(s2_l1c_ds.attrs, 'U')  # u_sun['U'] Reflectance_Conversion
tecta = util.interp_sza_s2(s2_l1c_ds, 60)
quant = get_value(s2_l1c_ds.attrs, "QUANTIFICATION_VALUE")
u_diff_temp = (get_value(s2_l1c_ds.attrs, 'DATASTRIP_SENSING_START') - datetime.datetime(2015, 6, 23, 10, 00)).days / 365.25 * 0.15

# Same value(s) for all S2A/S2B

# u_diff_cos = 0.4  # [%] from 0.13° diffuser planarity/micro as in (AIRBUS 2015). Assumed same for S2A/S2B.
# u_diff_k = 0.3  # [%] as a conservative residual (AIRBUS 2015). Assumed same for S2A/S2B.
# u_ADC =  # u_ADC = 0.5  # [DN](rectangular distribution, see combination)
# u_gamma =  # 0.4

# Running S2-RUT

B03 = s2_rut_algo.S2RutAlgo(a, e_sun, u_sun, tecta, quant, alpha, beta, 0.4,
                            0.3, u_diff_temp, 0.5, 0.4, 1)

# B03_err = B03.unc_calculation(s2_l1c_ds.B03.values, 2, 'Sentinel-2A')  # when roi is defined

B03_err = B03.unc_calculation(subset, 0, 'Sentinel-2A')

plt.imshow(B03_err)
plt.colorbar()
plt.show()

print(B03_err)

if __name__ == "__main__":
    pass
