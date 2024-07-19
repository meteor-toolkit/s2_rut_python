from typing import Optional, Union, List, Dict
import numpy as np
import s2_rut_algo
import datetime
from s2_l1_rad_conf import *
from eoio.utils.dict_tools import *
from eoio.processors import utils as util
import xarray as xr
import re  # RegEx checking strings for specific patterns

__author__ = [
    "Sam Hunt <sam.hunt@npl.co.uk>",
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
    "Maddie Stedman <maddie.stedman@npl.co.uk>",
]

__all__ = [
    "S2RUT"
]

meas_var_res = {
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

class S2RUT:
    DEFAULT_BAND_UNC_PARAMS = {
        'a': 0.0,
        'e_sun': 0.0,
        'u_sun': 1.0,
        'tecta': 0.0,
        'quant': 10000.0,
        'alpha': 0.0,
        'beta': 0.0,
        'u_diff_temp': 1.0,
    }

    DEFAULT_GLOBAL_UNC_PARAMS = {
        'u_diff_cos': 0.4,
        'u_diff_k': 0.3,
        'u_ADC': 0.5,
        'u_gamma': 0.4,
        'k': 1.0,
    }

    def run_rut(self,
                data_set: xr.Dataset,
                band_names: Union[list, str] = None,
                unc_contributions: Union[list, str] = None,
                unc_params: Union[list, str] = None,  # {'u_sun': , ..., 'k': ...}
                unc_correlations: Optional[Union[list, str]] = None,  # 'systematic', ['systematic', 'random']
                coverage_factor: int=1) -> Dict[str, Dict[str, float]]:
        """
        Run the Sentinel 2 radiometric uncertainty tool s2_rut

        :param data_set: satellite data_set product for which to calculate uncertainties
        :param band_names: definition of desired S2 bands, by default None
        :param unc_contributions: definition of the desired uncertainty contributions, by default None
        :param unc_params: definition of the desired uncertainty parameters, by default None
        :param unc_correlations:  definition of the desired uncertainty correlations, by default None
        :param coverage_factor: chosen coverage factor defining the confidence level of the uncertainty
        """

        # Validating input variables and giving errors if they are incorrect




        # select unc contributions to run

        # set unc params
        band_unc_params = self.get_band_unc_parameters(data_set, band_names)

        if band_unc_params is None:
            unc_params = {self.DEFAULT_GLOBAL_UNC_PARAMS, **self.DEFAULT_BAND_UNC_PARAMS}  # get from user input
        else:
            unc_params = {self.DEFAULT_GLOBAL_UNC_PARAMS, **band_unc_params}


        # make S2RUTALGO
        rut = s2_rut_algo.S2RutAlgo(**unc_params) # k=1, u_sun=1,
        unc = rut.unc_calculation(data_set, self.band_id, data_set.platform)

        return unc

    def return_unc_contributions(self, unc_correlations, unc_contributions):
        """
        Return dictionary of uncertainty contributions to include, by default return true for all.
        """
        unc_contributions = {"instrument_noise": True,
                            "oof_straylight_systematic": True,
                            "oof_straylight_random": True,
                            "crosstalk": True,
                            "adc_quantisation": True,
                            "ds_stability": True,
                            "gamma_knowledge": True,
                            "diffuser_absolute_knowledge": True,
                            "diffuser_temporal_knowledge": True,
                            "diffuser_straylight_residual": True,
                            "l1c_image_quantisation": True}

        #From s2_rut.py line 365. Filters which unc have been selected.
        # def get_unc_select(self):
        #     """ Parse the params """
        #
        #     # get param keys
        #     pkeys = self.params.keys()
        #
        #     # potential args
        #     args = [
        #         "ADC_quantisation",
        #         "Crosstalk",
        #         "Diffuser-absolute_knowledge",
        #         "Diffuser-cosine_effect",
        #         "Diffuser-straylight_residual",
        #         "Diffuser-temporal_knowledge",
        #         "DS_stability",
        #         "Gamma_knowledge",
        #         "Instrument_noise",
        #         "L1C_image_quantisation",
        #         "OOF_straylight-random",
        #         "OOF_straylight-systematic",
        #     ]
        #     self.args = args
        #
        #     # all contributions flag
        #     if "all_contribs" in pkeys:
        #         if self.params["all_contribs"]:
        #             return np.ones(12).astype(bool)
        #         else:
        #             return np.zeros(12).astype(bool)
        #
        #     else:
        #         # otherwise sort out the individual settings
        #         settings = []
        #         for arg in args:
        #             if arg in pkeys:
        #                 if self.params[arg]:
        #                     settings.append(True)
        #                 else:
        #                     settings.append(False)
        #             else:
        #                 # if setting not in there then set to False as default
        #                 settings.append(False)
        #         return np.array(settings)

        # According to Javie's paper, certain combinations of uncertainty are fully random, systematic, or structured
        if "random":
            unc_correlations = None
        elif "systematic":
            unc_correlations = None
        elif "structured":
            unc_correlations = None

        return unc_contributions

    def get_band_unc_parameters(self, data_set, band_names):

        # Defining the band names, accounting for both name types (with 0 and without)
        band_names = ['B0?1$',
                      'B0?2$',
                      'B0?3$',
                      'B0?4$',
                      'B0?5$',
                      'B0?6$',
                      'B0?7$',
                      'B0?8$',
                      'B0?8A$',
                      'B0?9$',
                      'B10',
                      'B11',
                      'B12']
        self.band_id = band_names.index(band_names) # indexing the band names, needed for extracting info from metadata

        # Extracting band uncertainty information from the provided data_set using eoio
        band_params = []
        for bands in band_names:
            band_params[bands] = {
                'a':  data_set[bands].PHYSICAL_GAINS,
                'e_sun': data_set[bands].SOLAR_IRRADIANCE['#text'],
                'u_sun': get_value(data_set.attrs, 'U'),
                'tecta': util.interp_sza_s2(data_set, meas_var_res[self.band_id[bands]]),
                'quant': get_value(data_set.attrs, "QUANTIFICATION_VALUE"),
                'alpha': data_set[bands].ALPHA,
                'beta': data_set[bands].BETA,
                'u_diff_temp': (get_value(data_set.attrs, 'DATASTRIP_SENSING_START') - datetime.datetime(2015, 6, 23, 10,
                                                                                                         00)).days / 365.25 * \
                          u_diff_temp_rate[data_set.platform][self.band_id[bands]],
            }
        return band_params

if __name__ == "__main__":
    pass
