"""
Code utilising eoio to read S2 L1 data product and apply s2-rut algorithm to calculate uncertainties of the product.

Data set selected is from the eoio documentation, read_examples_s2.py

Created 02/05/2024 by R.Ormane
"""
import os
import datetime

import matplotlib.pyplot as plt

import s2_rut_algo
from eoio.interface import (
    read,
    product_bounds,
    mid_lon_lat,
)
from eoio.utils.dict_tools import *
from eoio.processors import utils as util

__all__ = []

s2_l1c_filepaths = os.path.abspath(r"T:\ECO\EOServer\data\satellite\S2A_MSI\L1\PICS\Libya4"
                                   r"\S2A_MSIL1C_20190113T090331_N0207_R007_T34RGS_20190113T111955.SAFE")

s2_l1c_product_bounds = product_bounds(s2_l1c_filepaths)
s2_l1c_mid_lat_lon = mid_lon_lat(s2_l1c_filepaths)

s2_l1c_ds = read(
    s2_l1c_filepaths,
    subset_info={
        "meas_vars": ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"],
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

bandList = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']

product = 'Sentinel-2A'

# Data from Global Metadata (same for all bands in one product)

min_pxl = 1000
max_pxl = 1501
u_sun = get_value(s2_l1c_ds.attrs, 'U')  # u_sun['U'] Reflectance_Conversion
tecta = util.interp_sza_s2(s2_l1c_ds, 10)
quant = get_value(s2_l1c_ds.attrs, "QUANTIFICATION_VALUE")
u_diff_temp = (get_value(s2_l1c_ds.attrs, 'DATASTRIP_SENSING_START') - datetime.datetime(2015, 6, 23, 10, 00)).days / 365.25 * 0.04


for bandID, bands in enumerate(bandList):
    a = s2_l1c_ds[bands].PHYSICAL_GAINS  # physical gains
    e_sun = s2_l1c_ds[bands].SOLAR_IRRADIANCE['#text']  # solar irradiance
    alpha = s2_l1c_ds[bands].ALPHA  # noise
    beta = s2_l1c_ds[bands].BETA  # noise
    subset = s2_l1c_ds[bands].values[min_pxl:max_pxl, min_pxl:max_pxl]
    band = s2_rut_algo.S2RutAlgo(a, e_sun, u_sun, tecta[min_pxl:max_pxl, min_pxl:max_pxl], quant, alpha, beta, 0.4,
                                0.3, u_diff_temp, 0.5, 0.4, 1)
    bandUnc = band.unc_calculation(subset, bandID, 'Sentinel-2A')
    plt.imshow(bandUnc)
    plt.colorbar()
    plt.show()
    print(bandUnc.T)

if __name__ == "__main__":
    pass
