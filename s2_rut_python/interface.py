"""s2_rut_interface - Sentinel-2A and Sentinel-2B L1 uncertainty calculation class """

import datetime
import os
import sys
import warnings
from typing import Dict, List, Optional, Union

import obsarray  # type: ignore[import-untyped]
import xarray as xr
from processor_tools.processor_tools.utils.dict_tools import get_value  # type: ignore[import-untyped]

import s2_rut_python.sza_interp as util

THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
S2_RUT_DIRECTORY = os.path.join(THIS_DIRECTORY, "snap-rut", "src", "main", "python")

sys.path.insert(0, S2_RUT_DIRECTORY)
import s2_l1_rad_conf as conf  # type: ignore
import s2_rut_algo as srut  # type: ignore

__author__ = [
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
    "Sam Hunt <sam.hunt@npl.co.uk>",
    "Maddie Stedman <maddie.stedman@npl.co.uk>",
]

__all__ = ["S2RUT", "MyS2RUTAlgo"]

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
# Info from Gorroño, Javier, et al. "Providing uncertainty estimates of the Sentinel-2 top-of-atmosphere measurements for radiometric validation activities." European Journal of Remote Sensing 51.1 (2018): 650-666.
COMPONENTS = {
    "systematic": [
        "OOF_straylight-systematic",
        "Crosstalk",
        "Diffuser-straylight_residual",
        "Diffuser-absolute_knowledge",
        "Diffuser-temporal_knowledge",
        "DS_stability",  # correlated along dims spatial & temporal - structured along spectral dim
        "Diffuser-cosine_effect",
    ],
    "random": [
        "OOF_straylight-random",  # random in spectral & spatial, systematic in temporal
        "Instrument_noise",  # random in temporal dim. random in spectral & spatial dims "under the assumption that noise introduced by the post-amplification is not dominant"
        "ADC_quantisation",  #
        "Gamma_knowledge",  # random in spatial dim, fully correlated in temporal & spectral dims
        "L1C_image_quantisation",
    ],
    "structured": [],
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
    # Overwrite set parameter values from ESA's S2RUT tool with those extracted from product metadata
    def __init__(self, param_dict):
        super().__init__()

        for key, value in param_dict.items():
            if hasattr(self, key):
                setattr(self, key, value)


class S2RUT:
    """
    Class to run ESA's S2RUT tool for dataset to return dataset with uncertainty variables added using obsarray.
    """

    def __init__(self):
        # Define the band names and index the band names
        self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}
        # Initialise ESA's S2RUT class to access the default set parameters
        self.og_rut = srut.S2RutAlgo()

    def run(
        self,
        ds: xr.Dataset,
        band_names: Union[List[str], str, bool] = True,
        components: bool = True,
    ) -> xr.Dataset:
        """
        Run the Sentinel 2 radiometric uncertainty tool s2_rut

        :param ds: satellite dataset product for which to calculate uncertainties
        :param band_names: definition of desired S2 bands,
                           options: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12, by default None
        :param components: boolean - if True uncertainty provided as components systematic/random/structured. If False, total uncertainty provided.

        :return: input ds with uncertainty variables added using obsarray
        """
        unc_comp = self.return_unc_components()

        if band_names is True:
            band_names = [
                "B01",
                "B02",
                "B03",
                "B04",
                "B05",
                "B06",
                "B07",
                "B08",
                "B8A",
                "B09",
                "B10",
                "B11",
                "B12",
            ]

        for band in band_names:  # type: ignore[union-attr]
            if band not in ds:
                warnings.warn(
                    f"{band} data not in dataset so uncertainties not calculated for {band}."
                )
                continue
            band_unc_params = self.get_band_unc_parameters(ds, band)
            rut = MyS2RUTAlgo(band_unc_params)

            # Add parameters to ds for each band
            for comp in unc_comp.keys():
                rut.unc_select = list(unc_comp[comp].values())
                unc = rut.unc_calculation(
                    ds[band].values, self.band_id[band], ds.attrs.get("platform")
                )
                err_corr_def = [
                    {
                        "dim": [ds[band].dims],
                        "form": comp,
                        "params": [],
                        "units": ds[band].attrs.get("units"),
                    },
                ]
                ds.unc[band][f"u_{comp}_{band}"] = (
                    ds[band].dims,
                    unc,
                    {"err_corr": err_corr_def, "pdf_shape": "gaussian"},
                )

        return ds

    def return_unc_components(self):
        """
        Returns dictionary of uncertainty correlations, organising them as random, systematic, and structured.
        """

        comp_types_all = {}
        for unc_type in COMPONENTS:
            unc_comp = {key: False for key in U_CONTRIBUTIONS}
            for key in COMPONENTS[unc_type]:
                unc_comp[key] = True
                comp_dict = {unc_type: unc_comp}
                # add contribution dictionary to the correlation dictionary
                comp_types_all.update(comp_dict)

        return comp_types_all

    def get_band_unc_parameters(self, ds, band):
        """
        Extract band-specific uncertainty parameters from the provided data_set (eoio specific).
        """

        # Extract band uncertainty information (using eoio)
        band_params = {
            "a": get_value(ds[band].attrs["product_metadata"], "PHYSICAL_GAINS"),
            "e_sun": get_value(
                ds[band].attrs["product_metadata"]["SOLAR_IRRADIANCE"], "#text"
            ),
            "u_sun": get_value(ds.attrs, "U"),
            "tecta": util.interp_sza_s2(ds, MEAS_VAR_RES[band]),
            "quant": get_value(ds.attrs, "QUANTIFICATION_VALUE"),
            "alpha": get_value(ds[band].attrs["product_metadata"], "ALPHA"),
            "beta": get_value(ds[band].attrs["product_metadata"], "BETA"),
            "u_diff_cos": self.og_rut.u_diff_cos,
            "u_diff_k": self.og_rut.u_diff_k,
            "u_diff_temp": (
                get_value(ds.attrs, "DATASTRIP_SENSING_START")
                - TIME_INIT[ds.attrs["platform"]]
            ).days
            / 365.25
            * conf.u_diff_temp_rate[ds.platform][self.band_id[band]],
            "u_ADC": self.og_rut.u_ADC,
            "u_gamma": self.og_rut.u_gamma,
            "k": self.og_rut.k,
        }

        return band_params


if __name__ == "__main__":
    pass
