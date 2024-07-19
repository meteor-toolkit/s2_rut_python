import os
from eoio.interface import read
from s2_rut_interface import S2RUT

s2_l1c_filepaths = os.path.abspath(
    r"T:\ECO\EOServer\data\satellite\S2A_MSI\L1\PICS\Libya4\S2A_MSIL1C_20190113T090331_N0207_R007_T34RGS_20190113T111955.SAFE")

s2_l1c_ds = read(
    s2_l1c_filepaths,
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
)


# Example running S2_RUT interface

s2rut = S2RUT()

test = s2rut.run_rut(
    data_set=s2_l1c_ds,
    band_names=['B01', 'B09'],  # ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'],  # ['B01', "B05", "B09", "B10"],
    unc_contributions=[
                        "ADC_quantisation",
                        "Crosstalk",
                        "Diffuser-absolute_knowledge",
                        "Diffuser-cosine_effect",
                        "Diffuser-straylight_residual",
                        "Diffuser-temporal_knowledge",
                        "DS_stability",
                        "Gamma_knowledge",
                        "Instrument_noise",
                        "L1C_image_quantisation",
                        "OOF_straylight-random",
                        "OOF_straylight-systematic",
                        ],
    unc_correlations=["systematic", "random", "structured"],  # or "random" or "structured"
    coverage_factor=1.0,
)

print(test)
