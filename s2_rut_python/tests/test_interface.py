import unittest
import unittest.mock as mock

import numpy as np
import xarray as xr

from s2_rut_python.interface import COMPONENTS, MEAS_VAR_RES, U_CONTRIBUTIONS, S2RUTTool


class TestS2RUTTool(unittest.TestCase):
    def setUp(self):
        self.tool = S2RUTTool()
        self.ds = self._build_dataset()

    @staticmethod
    def _build_dataset() -> xr.Dataset:
        b01_data = np.array([[1.0, 0.0], [2.0, 3.0]], dtype=float)
        sza_data = np.array([[30.0, 30.0], [30.0, 30.0]], dtype=float)

        return xr.Dataset(
            data_vars={
                "B01": (
                    ["y_60m", "x_60m"],
                    b01_data,
                    {
                        "units": "1",
                        "long_name": "Reflectance in band B01",
                        "product_metadata": {
                            "solar_irradiance": 1800.0,
                            "noise_model_alpha": 0.4,
                            "noise_model_beta": 0.01,
                            "physical_gains": 1.2,
                            "radiometric_offset": 0.0,
                        },
                    },
                ),
                "solar_zenith_angle": (["y_60m", "x_60m"], sza_data),
            },
            coords={
                "y_60m": [0, 1],
                "x_60m": [0, 1],
            },
            attrs={
                "platform": "Sentinel-2A",
                "quantification_level": 10000,
                "reflectance_conversion_u": 0.05,
                "product_metadata": {
                    "platform": "Sentinel-2A",
                    "quantification_level": 10000,
                    "reflectance_conversion_u": 0.05,
                },
            },
        )

    def test___init__(self):
        self.assertEqual(self.tool.band_id["B01"], 0)
        self.assertEqual(len(self.tool.band_id), len(MEAS_VAR_RES))

    def test_set_contributor(self):
        class FakeRUT:
            unc_select_noise = False

        rut = FakeRUT()
        self.tool.set_contributor(rut, "noise", True)
        self.assertTrue(rut.unc_select_noise)

    def test_set_contributor_KeyError(self):
        class FakeRUT:
            pass

        with self.assertRaises(KeyError):
            self.tool.set_contributor(FakeRUT(), "not_a_contributor", True)

    def test__get_contributor_component_systematic(self):
        self.assertEqual(self.tool._get_contributor_component("u_sys"), "systematic")

    def test__get_contributor_component_random(self):
        self.assertEqual(self.tool._get_contributor_component("u_noise"), "random")

    def test__get_contributor_component_unknown(self):
        with self.assertRaises(ValueError):
            self.tool._get_contributor_component("u_unknown")

    def test__build_uncertainty_attrs_grouped(self):
        attrs = self.tool._build_uncertainty_attrs(self.ds, "B01", "grouped", "systematic")
        self.assertIn("Systematic", attrs["long_name"])
        self.assertEqual(attrs["units"], self.ds["B01"].attrs["units"])
        self.assertEqual(attrs["standard_name"], "uncertainty_systematic_B01")

    def test__build_uncertainty_attrs_per_contributor(self):
        attrs = self.tool._build_uncertainty_attrs(self.ds, "B01", "per_contributor", "u_noise")
        self.assertIn("Noise", attrs["long_name"])
        self.assertIn("u_noise", attrs["description"])
        self.assertEqual(attrs["units"], self.ds["B01"].attrs["units"])

    def test__compute_grouped_uncertainty_systematic(self):
        unc_contributors = {
            "u_sys": np.ones((2, 2)),
            "u_stray_sys": np.ones((2, 2)) * 2,
            "u_diff": np.ones((2, 2)) * 3,
        }
        result = self.tool._compute_grouped_uncertainty("systematic", unc_contributors)
        expected = np.ones((2, 2)) + np.sqrt((2 ** 2) + (3 ** 2))
        np.testing.assert_allclose(result, expected)

    def test__compute_grouped_uncertainty_random(self):
        unc_contributors = {
            "u_noise": np.ones((2, 2)) * 2,
            "u_adc": np.ones((2, 2)) * 3,
        }
        result = self.tool._compute_grouped_uncertainty("random", unc_contributors)
        expected = np.sqrt((2 ** 2) + (3 ** 2)) * np.ones((2, 2))
        np.testing.assert_allclose(result, expected)

    def test__configure_contributors_default(self):
        rut = mock.MagicMock()
        with mock.patch.object(self.tool, "set_contributor") as mocked:
            self.tool._configure_contributors(rut, None)

        self.assertEqual(mocked.call_count, len(U_CONTRIBUTIONS))
        for call in mocked.call_args_list:
            self.assertTrue(call.args[2])

    def test__configure_contributors_subset(self):
        rut = mock.MagicMock()
        subset = ["noise", "adc"]

        with mock.patch.object(self.tool, "set_contributor") as mocked:
            self.tool._configure_contributors(rut, subset)

        call_map = {c.args[1]: c.args[2] for c in mocked.call_args_list}
        self.assertTrue(call_map["noise"])
        self.assertTrue(call_map["adc"])
        self.assertFalse(call_map["ds"])

    def test__normalize_data_vars(self):
        self.assertEqual(self.tool._normalize_data_vars(True), list(MEAS_VAR_RES.keys()))
        self.assertEqual(self.tool._normalize_data_vars("B01"), ["B01"])
        self.assertEqual(self.tool._normalize_data_vars(["B01", "B09"]), ["B01", "B09"])

    def test__store_grouped_uncertainties(self):
        ds = self.ds.copy(deep=True)
        valid_mask = ds["B01"] != 0
        unc_contributors = {
            "u_sys": np.ones((2, 2)),
            "u_noise": np.ones((2, 2)) * 2,
            "u_stray_sys": np.ones((2, 2)) * 3,
            "u_diff": np.ones((2, 2)) * 4,
            "u_adc": np.ones((2, 2)) * 5,
        }

        names = self.tool._store_grouped_uncertainties(ds, "B01", unc_contributors, valid_mask)

        self.assertIn("u_systematic_B01", names)
        self.assertIn("u_random_B01", names)
        self.assertTrue(np.isnan(ds.unc["B01"]["u_random_B01"].value.values[0, 1]))
        self.assertEqual(ds.unc["B01"]["u_random_B01"].value.attrs["units"], ds["B01"].attrs["units"])

    def test__store_per_contributor_uncertainties(self):
        ds = self.ds.copy(deep=True)
        valid_mask = ds["B01"] != 0
        unc_contributors = {
            "u_noise": np.ones((2, 2)),
            "u_unknown": np.ones((2, 2)) * 2,
        }

        with mock.patch("s2_rut_python.interface.warnings.warn") as mocked_warn:
            names = self.tool._store_per_contributor_uncertainties(ds, "B01", unc_contributors, valid_mask)

        self.assertIn("u_noise_B01", names)
        self.assertIn("u_unknown_B01", names)
        self.assertTrue(np.isnan(ds.unc["B01"]["u_noise_B01"].value.values[0, 1]))
        mocked_warn.assert_called_once()

    def test__apply_zero_reflectance_mask(self):
        ds = self.ds.copy(deep=True)
        valid_mask = ds["B01"] != 0
        self.tool._apply_zero_reflectance_mask(ds, "B01", valid_mask)
        self.assertTrue(np.isnan(ds["B01"].values[0, 1]))
        self.assertEqual(ds["B01"].values[0, 0], 1.0)

    def test_return_sza_var(self):
        result = self.tool.return_sza_var(self.ds, "B01")
        self.assertEqual(result, "solar_zenith_angle")

    def test_return_sza_var_with_interp(self):
        ds = self.ds.copy(deep=True)
        ds["solar_zenith_angle"] = xr.DataArray(np.ones((1, 1)), dims=("a", "b"))
        ds["solar_zenith_angle_interp"] = xr.DataArray(np.ones((2, 2)), dims=("y_60m", "x_60m"))

        result = self.tool.return_sza_var(ds, "B01")
        self.assertEqual(result, "solar_zenith_angle_interp")

    def test_return_sza_var_missing(self):
        ds = self.ds.drop_vars("solar_zenith_angle")
        with self.assertRaises(KeyError):
            self.tool.return_sza_var(ds, "B01")

    def test_return_metadata(self):
        metadata = self.tool.return_metadata(self.ds, ["B01"])
        self.assertEqual(metadata["spacecraft"], "Sentinel-2A")
        self.assertIn("B01", metadata["Esun"])
        self.assertIn("B01", metadata["alpha"])

    def test_return_metadata_missing(self):
        ds = self.ds.copy(deep=True)
        ds.attrs.pop("platform", None)
        ds.attrs["product_metadata"].pop("platform", None)
        with self.assertRaises(KeyError):
            self.tool.return_metadata(ds, ["B01"])

    @mock.patch("s2_rut_python.interface.S2RUT_L1")
    def test_run_group_unc(self, mock_rut_cls):
        ds = self.ds.copy(deep=True)

        rut_instance = mock.MagicMock()
        for contributor in U_CONTRIBUTIONS:
            setattr(rut_instance, f"unc_select_{contributor}", True)

        rut_instance.unc_calculation_abs.return_value = (
            np.ones((2, 2)),
            {
                "u_sys": np.ones((2, 2)),
                "u_noise": np.ones((2, 2)) * 2,
                "u_stray_sys": np.ones((2, 2)) * 3,
                "u_diff": np.ones((2, 2)) * 4,
                "u_adc": np.ones((2, 2)) * 5,
            },
        )
        mock_rut_cls.return_value = rut_instance

        out = self.tool.run(ds=ds, data_vars=["B01"], group_unc=True)

        self.assertIn("u_random_B01", list(out.unc["B01"].keys()))
        self.assertIn("u_systematic_B01", list(out.unc["B01"].keys()))
        self.assertTrue(np.isnan(out["B01"].values[0, 1]))
        rut_instance.unc_calculation_abs.assert_called_once()

    @mock.patch("s2_rut_python.interface.S2RUT_L1")
    def test_run_per_contributor(self, mock_rut_cls):
        ds = self.ds.copy(deep=True)

        rut_instance = mock.MagicMock()
        for contributor in U_CONTRIBUTIONS:
            setattr(rut_instance, f"unc_select_{contributor}", True)

        rut_instance.unc_calculation_abs.return_value = (
            np.ones((2, 2)),
            {
                "u_noise": np.ones((2, 2)) * 2,
            },
        )
        mock_rut_cls.return_value = rut_instance

        out = self.tool.run(ds=ds, data_vars="B01", group_unc=False, subset_unc=["noise"])

        self.assertIn("u_noise_B01", list(out.unc["B01"].keys()))
        self.assertTrue(rut_instance.unc_select_noise)
        self.assertFalse(rut_instance.unc_select_ds)


if __name__ == "__main__":
    unittest.main()
