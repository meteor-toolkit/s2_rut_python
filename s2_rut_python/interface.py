"""s2_rut_interface - Sentinel-2A and Sentinel-2B L1 uncertainty calculation class """

import datetime
import os
import warnings
from typing import Any, Dict, List, Optional, Sequence, Union

import obsarray  # type: ignore[import-untyped]  # noqa: F401
import xarray as xr
import numpy as np
from processor_tools.utils.dict_tools import get_value

THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
import s2_rut_python._vendor  # noqa: F401, E402
from S2RUT import S2RUT_L1  # noqa: E402


__author__ = [
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
    "Sam Hunt <sam.hunt@npl.co.uk>",
    "Maddie Stedman <maddie.stedman@npl.co.uk>",
]

__all__ = ["S2RUTTool"]

INPUT_CONTRIBUTORS = os.path.abspath(
    os.path.join(
        os.path.dirname(THIS_DIRECTORY), "third-party", "Data", "unc_contributors.json"
    )
)
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
# Ref: Gorroño et al., Remote Sensing 9(2):178 (2017)
COMPONENTS = {
    "systematic": [
        "sys",
        "stray_sys",
        "stray_xtalk",
        "xtalk",
        "ds",
        "diff_temp",
        "diff",
    ],
    "random": [
        "noise",
        "adc",
        "gamma",
        "ref_quant",
        "geoloc",
    ],
}


U_CONTRIBUTIONS = [
    "noise",  # "Instrument_noise",
    "stray_sys",  # "OOF_straylight-systematic",
    "stray_rand",  # "OOF_straylight-random",
    "xtalk",  # "Crosstalk",
    "adc",  # "ADC_quantisation",
    "ds",  # "DS_stability",
    "gamma",  # "Gamma_knowledge",
    "diff_abs",  # "Diffuser-absolute_knowledge",
    "diff_temp",  # "Diffuser-temporal_knowledge",
    "diff_cos",  # "Diffuser-cosine_effect",
    "diff_sl",  # "Diffuser-straylight_residual",
    "ref_quant",  # "L1C_image_quantisation",
    "geoloc",
]

OUTPUT_CONTRIBUTIONS = {
    "noise",
    "sys",
    "stray_xtalk",
    "ds",
    "gamma",
    "adc",
    "diff",
    "ref_quant",
    "geoloc",
}


class S2RUTTool:
    """
    Class to run ESA's S2RUT tool for dataset to return dataset with uncertainty variables added using obsarray.
    """

    def __init__(self):
        # Define the band names and index the band names
        self.band_id = {band: index for index, band in enumerate(MEAS_VAR_RES.keys())}

    def set_contributor(self, rut, contributor, value):
        """
        Set the value of a contributor in the S2RUT class.

        :param rut: instance of ESA's S2RUT class.
        :param contributor: name of the contributor to set.
        :param value: boolean value to enable/disable the contributor.
        :raises KeyError: if contributor not found in S2RUT class.
        """
        attr_name = f"unc_select_{contributor}"
        if hasattr(rut, attr_name):
            setattr(rut, attr_name, value)
        else:
            raise KeyError(f"Contributor '{contributor}' not found in S2RUT class.")

    def _get_contributor_component(self, contributor: str) -> str:
        """
        Determine the correlation component (systematic or random) for a given contributor.

        :param contributor: contributor name (e.g., 'u_noise').
        :return: component type ('systematic' or 'random').
        :raises ValueError: if contributor cannot be categorized.
        """
        contrib_base = contributor.split("u_")[-1]
        if contrib_base in COMPONENTS["systematic"]:
            return "systematic"
        elif contrib_base in COMPONENTS["random"]:
            return "random"
        else:
            raise ValueError(
                f"Contributor {contrib_base} not categorized as systematic or random. "
                f"Defaulting to random."
            )

    def _build_uncertainty_attrs(
        self,
        ds: xr.Dataset,
        band: str,
        unc_type: str,
        component_or_contributor: str,
    ) -> Dict:
        """
        Build metadata attributes for uncertainty variables.

        :param ds: xarray Dataset containing band metadata.
        :param band: band name.
        :param unc_type: 'grouped' or 'per_contributor'.
        :param component_or_contributor: component name ('systematic'/'random') or contributor name.
        :return: dictionary of attributes suitable for obsarray uncertainty variables.
        """
        band_long_name = ds[band].attrs.get("long_name", f"Band {band}")
        band_units = ds[band].attrs.get("units", "")

        if unc_type == "grouped":
            long_name = f"{component_or_contributor.capitalize()} radiometric uncertainty for {band_long_name}"
            description = f"Radiometric uncertainty ({component_or_contributor}) from all {component_or_contributor} contributors for {band}"
        else:  # per_contributor
            contrib_name = (
                component_or_contributor.replace("u_", "").replace("_", " ").title()
            )
            long_name = f"{contrib_name} radiometric uncertainty for {band_long_name}"
            description = (
                f"Radiometric uncertainty ({component_or_contributor}) for {band}"
            )

        attrs = {
            "long_name": long_name,
            "description": description,
            "units": band_units,
            "standard_name": f"uncertainty_{component_or_contributor}_{band}",
        }
        return attrs

    def _compute_grouped_uncertainty(
        self, component: str, unc_contributors: Dict
    ) -> np.ndarray:
        """
        Compute grouped uncertainty for a correlation component (systematic or random).

        :param component: 'systematic' or 'random'.
        :param unc_contributors: dictionary of per-contributor uncertainties.
        :return: numpy array of grouped uncertainty values.
        """
        return np.sqrt(
            np.sum(
                [
                    unc_contributors[f"u_{contrib}"] ** 2
                    for contrib in COMPONENTS[component]
                    if f"u_{contrib}" in unc_contributors
                ],
                axis=0,
            )
        )

    def _configure_contributors(self, rut, subset_unc: Optional[Sequence[str]]) -> None:
        """
        Configure which uncertainty contributors to include in the RUT calculation.

        :param rut: S2RUT_L1 instance.
        :param subset_unc: List of contributors to include; None means all contributors.
        """
        if subset_unc is not None:
            for contributor in U_CONTRIBUTIONS:
                enabled = contributor in subset_unc
                self.set_contributor(rut, contributor, enabled)
        else:
            # Enable all contributors
            for contributor in U_CONTRIBUTIONS:
                self.set_contributor(rut, contributor, True)

    def _normalize_data_vars(
        self, data_vars: Optional[Union[List[str], str, bool]]
    ) -> List[str]:
        """
        Normalize data_vars parameter into a list of band names.

        :param data_vars: True (all bands), a string (single band), or a list of bands.
        :return: list of band names to process.
        """
        if data_vars is True:
            return list(MEAS_VAR_RES.keys())
        elif isinstance(data_vars, str):
            return [data_vars]
        else:
            return list(data_vars) if data_vars else []

    def _store_grouped_uncertainties(
        self,
        ds: xr.Dataset,
        band: str,
        unc_contributors: Dict,
        valid_mask: xr.DataArray,
    ) -> List[str]:
        """
        Store uncertainty variables grouped by correlation type (systematic/random).

        :param ds: xarray Dataset to update with uncertainty variables.
        :param band: band name being processed.
        :param unc_contributors: per-contributor uncertainty values.
        :param valid_mask: boolean mask for valid (non-zero) reflectance pixels.
        :return: list of uncertainty variable names created.
        """
        created_names = []
        for component in COMPONENTS:
            unc = self._compute_grouped_uncertainty(component, unc_contributors)
            err_corr_def = [
                {
                    "dim": [ds[band].dims],
                    "form": component,
                    "params": [],
                    "units": ds[band].attrs.get("units"),
                },
            ]
            masked_unc = np.where(valid_mask.values, unc, np.nan)
            unc_name = f"u_{component}_{band}"
            unc_attrs = self._build_uncertainty_attrs(ds, band, "grouped", component)
            ds.unc[band][unc_name] = (
                ds[band].dims,
                masked_unc,
                {
                    "err_corr": err_corr_def,
                    "pdf_shape": "gaussian",
                    **unc_attrs,
                },
            )
            created_names.append(unc_name)
        return created_names

    def _store_per_contributor_uncertainties(
        self,
        ds: xr.Dataset,
        band: str,
        unc_contributors: Dict,
        valid_mask: xr.DataArray,
    ) -> List[str]:
        """
        Store uncertainty variables per-contributor (not grouped).

        :param ds: xarray Dataset to update with uncertainty variables.
        :param band: band name being processed.
        :param unc_contributors: per-contributor uncertainty values.
        :param valid_mask: boolean mask for valid (non-zero) reflectance pixels.
        :return: list of uncertainty variable names created.
        """
        created_names = []
        for contributor, unc in unc_contributors.items():
            # Determine the correlation component type for this contributor.
            try:
                component = self._get_contributor_component(contributor)
            except ValueError as e:
                warnings.warn(str(e))
                component = "random"

            err_corr_def = [
                {
                    "dim": [ds[band].dims],
                    "form": component,
                    "params": [],
                    "units": ds[band].attrs.get("units"),
                },
            ]
            masked_unc = np.where(valid_mask.values, unc, np.nan)
            unc_name = f"{contributor}_{band}"
            unc_attrs = self._build_uncertainty_attrs(
                ds, band, "per_contributor", contributor
            )
            ds.unc[band][unc_name] = (
                ds[band].dims,
                masked_unc,
                {
                    "err_corr": err_corr_def,
                    "pdf_shape": "gaussian",
                    "contributor": contributor,
                    **unc_attrs,
                },
            )
            created_names.append(unc_name)
        return created_names

    def _apply_zero_reflectance_mask(
        self, ds: xr.Dataset, band: str, valid_mask: xr.DataArray
    ) -> None:
        """
        Apply zero-reflectance masking to the reflectance band.
        Note: Uncertainty variables are already masked during creation.

        :param ds: xarray Dataset to update.
        :param band: band name being processed.
        :param valid_mask: boolean mask for valid (non-zero) reflectance pixels.
        """
        ds[band] = ds[band].where(valid_mask)

    def run(
        self,
        ds: xr.Dataset,
        data_vars: Optional[Union[List[str], str, bool]] = True,
        group_unc: Optional[bool] = True,
        subset_unc: Optional[Sequence[str]] = None,
    ) -> xr.Dataset:
        """
        Run the Sentinel-2 radiometric uncertainty tool (S2-RUT) on a dataset.

        :param ds: xarray Dataset with reflectance bands and required metadata.
        :param data_vars: Band names to process. True (default) processes all available bands,
                         a string processes one band, or a list processes specific bands.
                         Options: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09, B10, B11, B12.
        :param group_unc: If True (default), uncertainties grouped by correlation type (systematic/random).
                         If False, uncertainties stored per-contributor.
        :param subset_unc: Optional list of contributors to include. If None, all contributors are used.
        :return: Input dataset with uncertainty variables added via obsarray interface.
        """
        # Initialize RUT engine and configure selected uncertainty contributors.
        rut = S2RUT_L1(INPUT_CONTRIBUTORS)
        self._configure_contributors(rut, subset_unc)

        # Normalize input and extract metadata.
        data_vars = self._normalize_data_vars(data_vars)
        metadata = self.return_metadata(ds, data_vars)

        # Process each band.
        for band in data_vars:  # type: ignore[union-attr]
            if band not in ds:
                warnings.warn(
                    f"{band} not in dataset; skipping uncertainty calculation."
                )
                continue

            # Validate solar zenith angle availability and compute valid pixel mask.
            solar_var = self.return_sza_var(ds, band)
            valid_mask = ds[band] != 0

            # Build sun_zenith list indexed by band index (for RUT tool compatibility).
            # sun_zenith[bandind] must correspond to solar zenith for band at MEAS_VAR_RES index.
            bandind = self.band_id[band]
            sun_zenith = [None] * len(MEAS_VAR_RES)
            sun_zenith[bandind] = ds[solar_var].values

            # Calculate uncertainties for this band.
            unc_total, unc_contributors = rut.unc_calculation_abs(
                ds[band].values,
                band,
                bandind,
                metadata,
                sun_zenith,
                do_contributor=True,
            )

            # Store uncertainties (grouped or per-contributor).
            if group_unc:
                self._store_grouped_uncertainties(
                    ds, band, unc_contributors, valid_mask
                )
            else:
                self._store_per_contributor_uncertainties(
                    ds, band, unc_contributors, valid_mask
                )

            # Apply zero-reflectance masking.
            self._apply_zero_reflectance_mask(ds, band, valid_mask)

        return ds

    def return_sza_var(self, ds: xr.Dataset, band: str) -> str:
        """
        Check that the solar zenith angle variable is present in the dataset and has a shape compatible with the shape of the band data, required for uncertainty calculation. If not, raise error.

        :param ds: dataset read in by eoio containing the required solar zenith angle variable for uncertainty calculation
        :param band: name of the band for which to check the solar zenith angle variable
        :return: name of the solar zenith angle variable to use for uncertainty calculation
        """

        # check shape of solar zenith angle variable is compatible with shape of band data, if not, raise error
        solar_var = "solar_zenith_angle"
        for solar_var in [
            "solar_zenith_angle",
            "solar_zenith_angle_interp",
            f"solar_zenith_angle_{get_value(ds[band].attrs, 'geometry_id')}",
        ]:
            if solar_var not in ds:
                continue
            elif ds[solar_var].shape == ds[band].shape:
                return solar_var
        raise KeyError(
            f"Required solar zenith angle variable not found with compatible shape for band {band} "
            f"with resolution {get_value(ds[band].attrs, 'geometry_id')}."
        )

    def return_metadata(
        self, ds: xr.Dataset, band_names: Sequence[str]
    ) -> Union[Dict, KeyError]:
        """
        Extract required metadata for uncertainty calculation from the dataset read in by eoio.

        :param ds: dataset read in by eoio containing the required metadata for uncertainty calculation
        :param band_names: list of band names for which to extract metadata
        :return: dictionary of metadata parameters required for uncertainty calculation
        """
        metadata: Dict[str, Any] = {
            "spacecraft": None,
            "quant": None,
            "A": [],
            "offset": [],
            "alpha": {},
            "beta": {},
            "Esun": [],
            "Usun": None,
            "refined": False,
        }  # this dictionary contains the relevant metadata parameters

        # read required metadata from ds.attrs and fill the metadata dictionary
        METADATA_MAP = {
            "spacecraft": "platform",
            "quant": "quantification_level",
            "Usun": "reflectance_conversion_u",
        }
        VAR_MTD_MAP = {
            "Esun": "solar_irradiance",
            "alpha": "noise_model_alpha",
            "beta": "noise_model_beta",
            "A": "physical_gains",
            "offset": "radiometric_offset",
        }

        for key, path in METADATA_MAP.items():
            if path in ds.attrs:
                metadata[key] = ds.attrs[path]
            elif path in ds.attrs["product_metadata"]:
                metadata[key] = ds.attrs["product_metadata"][path]
            else:
                raise KeyError(
                    f"Required metadata parameter '{key}' not found in dataset attributes at path '{path}'."
                )

        for key, path in VAR_MTD_MAP.items():
            try:
                metadata[key] = {
                    band: ds[band].attrs["product_metadata"][path]
                    for band in band_names
                    if band in ds.data_vars
                }
            except KeyError:
                raise KeyError(
                    f"Required metadata parameter '{key}' not found in dataset attributes at path '{path}'."
                )

        return metadata


if __name__ == "__main__":
    pass
