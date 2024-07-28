"""run_s2rut_example - Example usage script for interface functions"""

import os
from eoio.interface import read
from s2_rut_interface import S2RUT

__author__ = "Rasma Ormane <rasma.ormane@npl.co.uk>"

# get data directory
DATA_DIRECTORY = os.path.abspath(r"T:\ECO\EOServer\data\satellite")

# get satellite product filepaths for A and B
s2a_msi_filepath = os.path.join(DATA_DIRECTORY, "S2A_MSI", "L1")
s2b_msi_filepath = os.path.join(DATA_DIRECTORY, "S2B_MSI", "L1")

# Example product for S2A
s2a_product_path = os.path.join(s2a_msi_filepath,
                                "PICS",
                                "Libya4",
                                "S2A_MSIL1C_20190113T090331_N0207_R007_T34RGS_20190113T111955.SAFE")

# Example product for SBA
s2b_product_path = os.path.join(s2b_msi_filepath,
                                "PICS",
                                "Libya4",
                                "S2B_MSIL1C_20190108T090339_N0207_R007_T34RGS_20190108T110320.SAFE")

# read and convert a L1 S2A or S2B satellite product into reflectance using eoio
s2_l1c_ds = read(
    s2a_product_path,
    subset_info={
        "meas_vars": ['B01', 'B09'],  # ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'],
        "metadataLevel": True,
        "masks": None,
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
    read_params={
        "save_extracted": True,
    }
)


# example running S2_RUT interface, defining bands and uncertainty type of interest

s2rut = S2RUT()
test = s2rut.run_rut(
    data_set=s2_l1c_ds,
    band_names=['B01', 'B09'],  # ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
    unc_info='components',  # 'total'
    # unc_correlations=["systematic", "random", "structured"], # If none defined, only returns combined uncertainty
)

print(test)