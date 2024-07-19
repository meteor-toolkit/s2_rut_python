from typing import Union, Optional, List, Dict
import re
import xarray as xarray
import datetime
from eoio.utils.dict_tools import *
from eoio.processors import utils as util

import sys
import os

THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
S2_RUT_DIRECTORY = os.path.join(THIS_DIRECTORY, "snap-rut", "src", "main", "python")

sys.path.append(S2_RUT_DIRECTORY)

from s2_rut_algo import S2RutAlgo
from s2_l1_rad_conf import *

__author__ = [
    "Sam Hunt <sam.hunt@npl.co.uk>",
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
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

# band_patterns = ['B0?1$', 'B0?2$', 'B0?3$', 'B0?4$', 'B0?5$', 'B0?6$', 'B0?7$', 'B0?8$', 'B0?8A$', 'B0?9$', 'B10$',
# 'B11$', 'B12$']
# band_patterns = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']

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
                coverage_factor: Optional[float] = None,
                band_names: Union[List[str], str] = None,
                unc_contributions: Optional[Union[List[str], str]] = None,
                unc_params: Optional[Union[List[str], str]] = None,
                unc_correlations: Optional[Union[List[str], str]] = None,
                ) -> Dict[str, Dict[str, float]]:
        """
        Run the Sentinel 2 radiometric uncertainty tool s2_rut

        :param data_set: satellite data_set product for which to calculate uncertainties
        :param band_names: definition of desired S2 bands, by default None
        :param unc_contributions: definition of the desired uncertainty contributions, by default None
        :param unc_params: definition of the desired uncertainty parameters, by default None
        :param unc_correlations: definition of the desired uncertainty correlations, by default None
        :param coverage_factor: chosen coverage factor defining the confidence level of the uncertainty
        """

        # Validate band inputs
        if band_names is None:
            raise ValueError(f"No bands have been selected. Define one of the following {list(MEAS_VAR_RES.keys())}.")
        else:
            if isinstance(band_names, str):
                band_names = [band_names]
            if not all(isinstance(band, str) for band in band_names):
                raise ValueError("Band names must be a list of strings.")  # wrong input format

            # Make sure the input bands exist for the satellite/instrument

            valid_band_names = []
            for band in band_names:
                if any(re.match(pattern, band) for pattern in MEAS_VAR_RES.keys()):
                    valid_band_names.append(band)
                else:
                    raise ValueError(f"Error: Provided band {band} does not exist for the provided data. Valid "
                                     f"bands are {list(MEAS_VAR_RES.keys())}")

        unc_info = {}
        unc_contr = self.return_unc_contributions(unc_contributions)
        unc_corr = self.return_unc_correlations(unc_correlations)

        # Both: uncertainty contributions and correlations are defined
        if unc_contributions is not None and unc_correlations is not None:
            unc_info.update(unc_corr)
            unc_info.update(unc_contr)
        # Only contributions are defined
        elif unc_contributions is not None:
            unc_info.update(unc_contr)
        # Only correlations are defined
        elif unc_correlations is not None:
            unc_info.update(unc_corr)
        # No contributions or correlations are defined
        else:
            unc_info.update(unc_contr)
            print(f"Warning: No uncertainty correlation type or contributions specified. By default all selected: '{U_CONTRIBUTIONS}'.")

        # If coverage factor is not 1.0, 2.0 or 3.0, set it to 1.0 by default
        if coverage_factor not in (1.0, 2.0, 3.0):
            coverage_factor = self.DEFAULT_GLOBAL_UNC_PARAMS['k']
            ValueError("Coverage_factor must equal to 1.0, 2.0 or 3.0. By default k has been set to 1.0.")

        # Obtain information about user selected contributions and/or correlations

        # Define a dictionary to store xarrays for each band
        datasets = {}
        for band in band_names:
            datasets[band] = xarray.Dataset()
            # 
            datasets[band]['uncertainty_variables'] = xarray.DataArray(data=unc_info)
            # Retrieve band uncertainty parameters from the provided data set and k for each band
            band_unc_params = self.get_band_unc_parameters(data_set, band, coverage_factor)

            # If no parameters are defined, use default. Else combine global and given parameters
            if unc_params is None:
                unc_params = {**self.DEFAULT_GLOBAL_UNC_PARAMS, **self.DEFAULT_BAND_UNC_PARAMS, **band_unc_params}
            else:
                unc_params = {**self.DEFAULT_GLOBAL_UNC_PARAMS, **band_unc_params, **unc_params}

            # Add parameters to xarray for each band then clear the unc_params for next band
            datasets[band]['parameters'] = xarray.DataArray(data=unc_params)
            unc_params = None

            # Initialize and run the uncertainty calculation algorithm for each of the bands
            for unc_type in unc_info.keys():
                rut = S2RutAlgo(**datasets[band].parameters.item())
                # overwrite Javie's unc_select with the user selected uncertainty parameter
                rut.unc_select = list(datasets[band].uncertainty_variables.item()[unc_type].values())
                unc = rut.unc_calculation(data_set[band].values, self.band_id[band], data_set.platform)
                datasets[band][unc_type] = xarray.DataArray(
                    data=unc,
                    dims=("x", "y"),
                )

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

    def return_unc_contributions(self, unc_contributions):
        """
        Return dictionary of uncertainty contributions to include, by default return true for all.
        """

        default_unc_contributions = {}
        # If no uncertainty contributions defined, all are selected
        if unc_contributions is None:
            for unc in U_CONTRIBUTIONS:
                default_unc_contributions[unc] = True
        # Filter which uncertainties have been selected
        elif unc_contributions is not None:
            for unc in U_CONTRIBUTIONS:
                if unc in unc_contributions:
                    default_unc_contributions[unc] = True
                else:
                    default_unc_contributions[unc] = False

        # Formats the contribution dict same as for correlations
        default_unc_contributions_dict = {'combined': default_unc_contributions}

        return default_unc_contributions_dict

    def get_band_unc_parameters(self, data_set, band, coverage_factor):
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
            'u_diff_temp': (get_value(data_set.attrs, 'DATASTRIP_SENSING_START')
                            - TIME_INIT[data_set.platform]).days / 365.25
                           * u_diff_temp_rate[data_set.platform][self.band_id[band]],
            'k': coverage_factor
        }

        return band_params


if __name__ == "__main__":
    pass
