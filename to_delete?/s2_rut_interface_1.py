"""s2_rut_interface - Sentinel-2A and Sentinel-2B L1 uncertainty calculation class """

from typing import Union, Optional, List, Dict
import re
import xarray as xarray
import datetime
from eoio.utils.dict_tools import *
from eoio.processors import utils as util
import sys
import os
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
# s2 = rut.S2RutAlgo()
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

# class MyS2RUTAlgo(srut.S2RutAlgo):
#     # Overwrite Javie's set parameter values with those extracted from metadata
#     def __init__(self, a, e_sun, u_sun, tecta, quant, alpha, beta, u_diff_cos, u_diff_k, u_diff_temp, u_ADC, u_gamma, k):
#         self.a = a
#         self.e_sun = e_sun
#         self.u_sun = u_sun
#         self.tecta = tecta
#         self.quant = quant
#         self.alpha = alpha
#         self.beta = beta
#         self.u_diff_cos = u_diff_cos
#         self.u_diff_k = u_diff_k
#         self.u_diff_temp = u_diff_temp
#         self.u_ADC = u_ADC
#         self.u_gamma = u_gamma
#         self.k = k
#         self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}
#


class MyS2RUTAlgo(srut.S2RutAlgo):
    # Overwrite Javie's set parameter values with those extracted from metadata
    def __init__(self, **kwargs):
        super().__init__()

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

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
                data_set: xarray.Dataset,
                band_names: Union[List[str], str] = None,
                unc_info: str = None,
                unc_correlations: Optional[Union[List[str], str]] = None,
                ) -> Dict[str, Dict[str, float]]:
        """
        Run the Sentinel 2 radiometric uncertainty tool s2_rut

        :param data_set: satellite data_set product for which to calculate uncertainties
        :param band_names: definition of desired S2 bands,
                           options: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12, by default None
        :param unc_info: chosen output uncertainty information
                         options: full, components, total
        :param unc_correlations: definition of the desired uncertainty correlations, by default None
        """

        unc_info_dict = {}
        unc_corr = self.return_unc_correlations(unc_correlations)

        # Correlations are defined
        if unc_correlations is not None:
            unc_info_dict.update(unc_corr)
        # No contributions or correlations are defined
        else:
            print(f"Warning: No uncertainty correlation typed. By default all selected: '{U_CONTRIBUTIONS}'.")
            unc_info_dict.update(self.corr_types.keys())

        # Define desired uncertainty information: full, components, total
        unc_info_options = ['components', 'total']
        # Desired uncertainty output has not been defined
        if unc_info is None:
            print(
                f"Warning: No uncertainty options specified. By default total uncertainty returned.")
            unc_info = unc_info_options
        for info in unc_info:
            if info not in unc_info_options:
                print(f"Warning: Invalid uncertainty options specified. Options are '{unc_info_options}'. By default total uncertainty returned.")
                unc_info = unc_info_options

        datasets = {}
        for band in band_names:
            datasets[band] = xarray.Dataset()
            band_unc_params = self.get_band_unc_parameters(data_set, band)
            if "full" in unc_info:
                # Add parameters to xarray for each band then clear the unc_params for next band
                datasets[band]['full'] = xarray.DataArray(data=band_unc_params)

            elif "total" in unc_info:
                unc_contr = self.return_unc_contributions()
                rut = MyS2RUTAlgo(band_unc_params)
                rut.unc_select = list(unc_contr.values())
                unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
                datasets[band]['total'] = xarray.DataArray(
                    data=unc,
                    dims=("x", "y"),
                )

            # if "components" in unc_info:
            #     unc_corr = self.return_unc_correlations(unc_correlations)
            # # Initialize and run the uncertainty calculation algorithm for each of the bands
            # for unc_type in unc_info_dict.keys():
            #     rut = MyS2RUTAlgo(**datasets[band].parameters.item())
            #     # overwrite Javie's unc_select with the user selected uncertainty parameter
            #     rut.unc_select = list(datasets[band].uncertainty_variables.item()[unc_type].values())
            #     unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
            #     datasets[band][unc_type] = xarray.DataArray(
            #         data=unc,
            #         dims=("x", "y"),
            #     )
            #
            # if "total" in unc_info:
            #     unc_contr = self.return_unc_contributions(unc_contributions)
            #             #     rut = MyS2RUTAlgo(**datasets[band].parameters.item())
            #             #     rut.unc_select = list(unc_contr.values())
            #             #     unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
            #             #     datasets[band]['total'] = xarray.DataArray(
            #             #         data=unc,
            #             #         dims=("x", "y"),
            #             #     )

        return datasets

    def return_unc_correlations(self, unc_correlations):
        """
        Return dictionary of uncertainty contributions to include, by default return true for all.
        """

        # Categorise correlations for all possible uncertainty contributions
        corr_types = {"systematic": ["OOF_straylight-systematic",
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

        corr_types_all = {}
        # If uncertainty correlations are defined, output the required contributions
        if unc_correlations is not None:
            for unc_type in unc_correlations:
                if unc_type in corr_types:
                    selected_unc_corr = {key: False for key in U_CONTRIBUTIONS}
                    for key in corr_types[unc_type]:
                        selected_unc_corr[key] = True
                    corr_dict = {unc_type: selected_unc_corr}
                    # add contribution dictionary to the correlation dictionary
                    corr_types_all.update(corr_dict)
                    # clear the correlation dictionary for the next correlation type
                    corr_dict.clear()
                else:
                    # If the provided correlation type is not recognised/misspelled give warning
                    print(
                        f"Warning: Unrecognised correlation '{unc_type}' specified. Available types are '{list(corr_types.keys())}'.")

        return corr_types_all

    def return_unc_contributions(self):
        """
        Return dictionary of uncertainty contributions to include, by default return true for all.
        """

        default_unc_contributions = {}
        for unc in U_CONTRIBUTIONS:
            default_unc_contributions[unc] = True

        return default_unc_contributions


    def get_band_unc_parameters(self, data_set, band):
        """
        Extract band-specific uncertainty parameters from the provided data_set (eoio specific).
        """
        # Define the band names and index the band names
        self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}

        # Extract band uncertainty information
        band_params = {
            'a': data_set[band].PHYSICAL_GAINS,
            'e_sun': data_set[band].SOLAR_IRRADIANCE['#text'],
            'u_sun': get_value(data_set.attrs, 'U'),
            'tecta': util.interp_sza_s2(data_set, MEAS_VAR_RES[band]),
            'quant': get_value(data_set.attrs, "QUANTIFICATION_VALUE"),
            'alpha': data_set[band].ALPHA,
            'beta': data_set[band].BETA,
            'u_diff_cos': self.DEFAULT_GLOBAL_UNC_PARAMS['u_diff_cos'],
            'u_diff_k': self.DEFAULT_GLOBAL_UNC_PARAMS['u_diff_k'],
            'u_diff_temp': (get_value(data_set.attrs, 'DATASTRIP_SENSING_START')
                            - TIME_INIT[data_set.platform]).days / 365.25
                           * conf.u_diff_temp_rate[data_set.platform][self.band_id[band]],
            'u_ADC': self.DEFAULT_GLOBAL_UNC_PARAMS['u_ADC'],
            'u_gamma': self.DEFAULT_GLOBAL_UNC_PARAMS['u_gamma'],
            'k': self.DEFAULT_GLOBAL_UNC_PARAMS['k']
        }

        return band_params


if __name__ == "__main__":
    pass
