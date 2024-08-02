"""s2_rut_interface - Sentinel-2A and Sentinel-2B L1 uncertainty calculation class """

from typing import Union, Optional, List, Dict
import xarray as xarray
import datetime
from eoio.utils.dict_tools import *
from eoio.processors import utils as util
import sys
import os
import numpy as np
import obsarray

THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
S2_RUT_DIRECTORY = os.path.join(THIS_DIRECTORY, "snap-rut", "src", "main", "python")

sys.path.insert(0, S2_RUT_DIRECTORY)
import s2_rut_algo as srut
import s2_l1_rad_conf as conf

__author__ = [
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
    "Sam Hunt <sam.hunt@npl.co.uk>",
    "Maddie Stedman <maddie.stedman@npl.co.uk>",
]

__all__ = [
    "S2RUT"
]

MEAS_VAR_RES = {
    "B01": 60,
    "B02": 10,
    "B03": 10,
    "B04": 10,
    "B05": 20,
    "B06": 20,
    "B07": 20,
    "B08": 10,
    "B8A": 20,
    "B09": 60,
    "B10": 60,
    "B11": 20,
    "B12": 20,
}

TIME_INIT = {
    "Sentinel-2A": datetime.datetime(2015, 6, 23, 10, 00),
    "Sentinel-2B": datetime.datetime(2017, 3, 7, 10, 00),
}

# Categorise correlations for all possible uncertainty contributions
CORR_TYPES = {"systematic": ["OOF_straylight-systematic",
                             "Crosstalk",
                             "Diffuser-straylight_residual",
                             "Diffuser-temporal_knowledge"
                             ],
              "random": ["OOF_straylight-random"],
              "structured": ["Instrument_noise",
                             "ADC_quantisation",
                             "DS_stability",
                             "Gamma_knowledge",
                             "Diffuser-absolute_knowledge",
                             "Diffuser-cosine_effect",
                             "L1C_image_quantisation",
                             ],
              }


U_CONTRIBUTIONS = [
    "Instrument_noise",
    "OOF_straylight-systematic",
    "OOF_straylight-random",
    "Crosstalk",
    "ADC_quantisation",
    "DS_stability",
    "Gamma_knowledge",
    "Diffuser-absolute_knowledge",
    "Diffuser-temporal_knowledge",
    "Diffuser-cosine_effect",
    "Diffuser-straylight_residual",
    "L1C_image_quantisation",
]


class MyS2RUTAlgo(srut.S2RutAlgo):
    # Overwrite Javie's set parameter values with those extracted from metadata
    def __init__(self, param_dict):
        super().__init__()

        for key, value in param_dict.items():
            if hasattr(self, key):
                setattr(self, key, value)


class S2RUT:
    def __init__(self):
        # Define the band names and index the band names
        self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}
        # Initialise Javie's original class to access the default set parameters
        # (u_diff_cos, u_diff_k, u_ADC, u_gamma, u_k),
        self.og_rut = srut.S2RutAlgo()

    def run_rut(self,
                data_set: xarray.Dataset,
                band_names: Union[List[str], str] = None,
                unc_info: Optional[str] = None,
                ) -> Dict[str, Dict[str, float]]:
        """
        Run the Sentinel 2 radiometric uncertainty tool s2_rut

        :param data_set: satellite data_set product for which to calculate uncertainties
        :param band_names: definition of desired S2 bands,
                           options: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12, by default None
        :param unc_info: chosen output uncertainty information
                         options: components, total
        """
        unc_corr = self.return_unc_correlations()
        result = {}

        for band in band_names:
            y_dim = data_set[band].coords.dims[0]
            x_dim = data_set[band].coords.dims[1]
            y_coord = data_set[band].coords[y_dim].values
            x_coord = data_set[band].coords[x_dim].values

            band_unc_params = self.get_band_unc_parameters(data_set, band)
            rut = MyS2RUTAlgo(band_unc_params)

            dataset = xarray.Dataset(coords={y_dim: y_coord, x_dim: x_coord})

            # If the user wants individual uncertainty components
            if unc_info == 'components':
                # Add parameters to xarray for each band
                for corr in unc_corr.keys():
                    rut.unc_select = list(unc_corr[corr].values())
                    unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
                    dataset[corr] = ([y_dim, x_dim], unc)

            # If the user wants total uncertainty
            elif unc_info == 'total':
                rut.unc_select = self.og_rut.unc_select
                unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
                dataset['total'] = ([y_dim, x_dim], unc)
            else:
                raise ValueError(f"Warning: Invalid input: {unc_info}. Options are 'total' and 'components'.")

            result[band] = dataset

        return result

    def return_unc_correlations(self):
        """
        Return dictionary of uncertainty contributions to include, by default return true for all.
        """

        corr_types_all = {}
        for unc_type in CORR_TYPES:
            unc_corr = {key: False for key in U_CONTRIBUTIONS}
            for key in CORR_TYPES[unc_type]:
                unc_corr[key] = True
                corr_dict = {unc_type: unc_corr}
                # add contribution dictionary to the correlation dictionary
                corr_types_all.update(corr_dict)

        return corr_types_all

    def get_band_unc_parameters(self, data_set, band):
        """
        Extract band-specific uncertainty parameters from the provided data_set (eoio specific).
        """
        # self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}

        # Extract band uncertainty information (using eoio)
        band_params = {
            'a': data_set[band].PHYSICAL_GAINS,
            'e_sun': data_set[band].SOLAR_IRRADIANCE['#text'],
            'u_sun': get_value(data_set.attrs, 'U'),
            'tecta': util.interp_sza_s2(data_set, MEAS_VAR_RES[band]),
            'quant': get_value(data_set.attrs, "QUANTIFICATION_VALUE"),
            'alpha': data_set[band].ALPHA,
            'beta': data_set[band].BETA,
            'u_diff_cos': self.og_rut.u_diff_cos,
            'u_diff_k': self.og_rut.u_diff_k,
            'u_diff_temp': (get_value(data_set.attrs, 'DATASTRIP_SENSING_START')
                            - TIME_INIT[data_set.platform]).days / 365.25
                           * conf.u_diff_temp_rate[data_set.platform][self.band_id[band]],
            'u_ADC': self.og_rut.u_ADC,
            'u_gamma': self.og_rut.u_gamma,
            'k': self.og_rut.k
        }

        return band_params


if __name__ == "__main__":
    pass
