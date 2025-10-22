"""run_s2rut_example - Example usage script for interface functions"""

import os
import socket
import xarray as xr
# from eoio.interface import read
# import eoio.utils.read_utils as utils
# from processor_tools.utils.dict_tools import get_value
from s2_rut_python.interface import S2RUT
import datetime

__author__ = "Rasma Ormane <rasma.ormane@npl.co.uk>"

# get data directory
if socket.gethostname() == 'lyon.npl.co.uk' or socket.gethostname() == 'eoserver.npl.co.uk':
    DATA_DIRECTORY = os.path.abspath(r"/mnt/t/data/satellite")
else:
    DATA_DIRECTORY = os.path.abspath(r"T:\ECO\EOServer\data\satellite")

# get satellite product filepaths for A and B
s2a_msi_filepath = os.path.join(DATA_DIRECTORY, "S2A_MSI", "L1")
s2b_msi_filepath = os.path.join(DATA_DIRECTORY, "S2B_MSI", "L1")

# Example product for S2A
s2a_product_path = os.path.join(s2a_msi_filepath,
                                "PICS",
                                "Libya4",
                                "S2A_MSIL1C_20190113T090331_N0207_R007_T34RGS_20190113T111955.SAFE")

# Example product for S2B
s2b_product_path = os.path.join(s2b_msi_filepath,
                                "PICS",
                                "Libya4",
                                "S2B_MSIL1C_20190108T090339_N0207_R007_T34RGS_20190108T110320.SAFE")
#
# # # read and convert a L1 S2A or S2B satellite product into reflectance using eoio
# s2_l1c_ds = read(
#     s2a_product_path,
#     subset_info={
#         "meas_vars": ['B01', 'B09'],  # ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'],
#         "metadataLevel": 'partial',
#         "masks": None,
#         "aux_data": [
#             "observation_geometry",
#         ],
#     },
#     process_params={
#         "convert": ["reflectance"],
#         # "angles": {"angle_type": 'both'},
#     },
#     read_params={
#         "save_extracted": True,
#     }
# )
#
# for var in s2_l1c_ds:
#     print(var, s2_l1c_ds[var].attrs)
#     s2_l1c_ds[var].attrs = {}
#
# U = get_value(s2_l1c_ds.attrs, "U")
# DATASTRIP_SENSING_START = get_value(s2_l1c_ds.attrs, "DATASTRIP_SENSING_START")
# # s2_l1c_ds = utils.clean_and_convert_attrs(s2_l1c_ds)
#
# print('U', U)
# print('DATASTRIP_SENSING_START', DATASTRIP_SENSING_START)
# s2_l1c_ds.attrs = {}
#
# s2_l1c_ds.to_netcdf(f'{os.path.basename(s2a_product_path).split(".")[0]}.nc')
# kill
s2_l1c_ds = xr.open_dataset(f'{os.path.basename(s2a_product_path).split(".")[0]}.nc')

s2_l1c_ds.attrs = {
    'U':1.03412680819686,
    'DATASTRIP_SENSING_START': datetime.datetime(2019, 1, 13, 9, 6, 47),
    'platform': 'Sentinel-2A',
    'QUANTIFICATION_VALUE': 10000
}
s2_l1c_ds['B01'].attrs = {'scale_factor': 1.0, 'add_offset': 0.0, 'units': '', 'standard_name': 'toa_bidirectional_reflectance_B01', 'long_name': 'Reflectance in band B01', 'ancillary_vars': [], 'measurand': 'reflectance', 'geometry': '60m', 'resolution': 60.0, 'product_metadata': {'ABSOLUTE_CALIBRATION_ACCURACY': 5, 'ALPHA': 0.56, 'BETA': 0.00054, 'COMPRESSION_RATE': 2.655, 'CROSS_BAND_CALIBRATION_ACCURACY': 3, 'INTEGRATION_TIME': {'@unit': 'ms', '#text': 7.4473767}, 'MULTI_TEMPORAL_CALIBRATION_ACCURACY': 1, 'PHYSICAL_GAINS': 4.11350657, 'RESOLUTION': 60, 'SOLAR_IRRADIANCE': {'@unit': 'W/m²/µm', '#text': 1884.69}, 'rsr_units': '1e-6 m', 'response': [0.0, 0.0, 0.0017757417, 0.0040730606, 0.0036261429, 0.0035151988, 0.0057291626, 0.0037802919, 0.002636732, 0.0012621128, 0.001987583, 0.0013689131, 0.001250444, 0.0004634544, 0.000814293, 0.0013764316, 0.0014850859, 0.0018237347, 0.0016268175, 0.004392062, 0.029008098, 0.11874593, 0.32387507, 0.57281923, 0.7147275, 0.7619678, 0.78929704, 0.80862385, 0.81089383, 0.8241988, 0.8541581, 0.8707909, 0.887311, 0.92619926, 0.9822815, 1.0, 0.9752382, 0.9359634, 0.8899715, 0.8502105, 0.8256945, 0.7839024, 0.6141742, 0.3300711, 0.124108315, 0.04365694, 0.014749595, 0.0, 0.0], 'response_wavelengths': [0.41, 0.411, 0.412, 0.413, 0.414, 0.415, 0.416, 0.417, 0.418, 0.419, 0.42, 0.421, 0.422, 0.423, 0.424, 0.425, 0.426, 0.427, 0.428, 0.429, 0.43, 0.431, 0.432, 0.433, 0.434, 0.435, 0.436, 0.437, 0.438, 0.439, 0.44, 0.441, 0.442, 0.443, 0.444, 0.445, 0.446, 0.447, 0.448, 0.449, 0.45, 0.451, 0.452, 0.453, 0.454, 0.455, 0.456, 0.457, 0.458]}, 'band_central_wavelength': 0.4426950483545459}
s2_l1c_ds['B09'].attrs = {'scale_factor': 1.0, 'add_offset': 0.0, 'units': '', 'standard_name': 'toa_bidirectional_reflectance_B09', 'long_name': 'Reflectance in band B09', 'ancillary_vars': [], 'measurand': 'reflectance', 'geometry': '60m', 'resolution': 60.0, 'product_metadata': {'ABSOLUTE_CALIBRATION_ACCURACY': 5, 'ALPHA': 0.563, 'BETA': 0.10258, 'COMPRESSION_RATE': 2.655, 'CROSS_BAND_CALIBRATION_ACCURACY': 3, 'INTEGRATION_TIME': {'@unit': 'ms', '#text': 7.593408}, 'MULTI_TEMPORAL_CALIBRATION_ACCURACY': 1, 'PHYSICAL_GAINS': 8.50073103, 'RESOLUTION': 60, 'SOLAR_IRRADIANCE': {'@unit': 'W/m²/µm', '#text': 812.92}, 'rsr_units': '1e-6 m', 'response': [0.0, 0.0, 0.016629528, 0.061118573, 0.17407094, 0.38946456, 0.6645915, 0.87454116, 0.93695986, 0.96751016, 0.9893391, 0.9951269, 1.0, 0.9784576, 0.9806912, 0.9922335, 0.9879838, 0.99428314, 0.9834804, 0.97820014, 0.9502337, 0.952996, 0.9224031, 0.8557383, 0.70970225, 0.46429542, 0.21538426, 0.065341204, 0.016255962, 0.0, 0.0], 'response_wavelengths': [0.93, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.94, 0.941, 0.942, 0.943, 0.944, 0.945, 0.946, 0.947, 0.948, 0.949, 0.9499999, 0.951, 0.952, 0.953, 0.95400006, 0.95500004, 0.95600003, 0.957, 0.958, 0.959, 0.96]}, 'band_central_wavelength': 0.9450544655757902}


# example running S2_RUT interface, defining bands and uncertainty type of interest
s2rut = S2RUT()
ds = s2rut.run(
    ds=s2_l1c_ds,
    band_names=['B01', 'B09'],  # ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
)

print(ds)


if __name__ == '__main__':
    pass