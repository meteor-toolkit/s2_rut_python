"""
Example script to demonstrate use of sat2xr class to read in Sentinel-2 data
"""

from sat2xr import *

# filepath to S2 .SAFE folder
base = r'C:\Users\mg13\OneDrive - National Physical ' \
       r'Laboratory\Documents\Sep2022\getting_to_grips\testdata' \
       r'\S2B_MSIL1C_20220910T152639_N0400_R025_T18QYF_20220910T195115.SAFE'

# desired bands
bands = ['B01', 'B02', 'B8A']

# initialise object
obj = S2toxr(base, bands=bands)

# read in satellite metadata and image data into xarray
obj.read_meta()
obj.read_img()

# print object class documentation
print(obj.__doc__)

# to obtain quality and cloud mask variables run obj.get_masks()

# access populated xarray (obj.sat_ds) by setting as variable
xds = obj.sat_ds

print('Object is happy :)')