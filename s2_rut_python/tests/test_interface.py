import datetime

# import os, sys
import unittest
import unittest.mock as mock
from unittest.mock import MagicMock

import numpy as np
import xarray as xr

from s2_rut_python.interface import S2RUT

# THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
# S2_RUT_DIRECTORY = os.path.join(THIS_DIRECTORY, "snap-rut", "src", "main", "python")

# sys.path.insert(0, S2_RUT_DIRECTORY)
# import s2_rut_algo
# import s2_l1_rad_conf

__author__ = [
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
]


class TestS2RUT(unittest.TestCase):
    def setUp(self):
        # Create a mock dataset
        self.mock_dataset = xr.Dataset(
            data_vars=dict(
                B01=(
                    ["y_60m", "x_60m"],
                    np.ones((5, 5)),
                    {
                        "PHYSICAL_GAINS": 1,
                        "SOLAR_IRRADIANCE": {"#text": 1},
                        "ALPHA": 1,
                        "BETA": 1,
                    },
                )
            ),
            coords=dict(
                y=(["y_60m"], np.ones(5)),
                x=(["x_60m"], np.ones(5)),
            ),
            attrs={
                "platform": "Sentinel-2A",
                "U": 0.05,
                "QUANTIFICATION_VALUE": 10000,
                "DATASTRIP_SENSING_START": datetime.datetime(
                    2015, 6, 1, 10, 00
                ),  # np.datetime64("2022-06-01T10:00:00")
            },
        )

        self.s2_rut = S2RUT()

    @mock.patch("processor_tools.utils.dict_tools.get_value")
    @mock.patch("s2_rut_python.interface.util.interp_sza_s2")
    @mock.patch("s2_rut_algo.S2RutAlgo")
    def test_get_band_unc_parameters(
        self, mock_s2_rut_algo, mock_interp_sza_s2, mock_get_value
    ):
        # # Set up the mock return values for get_value and interp_sza_s2
        mock_get_value.side_effect = lambda attrs, key: attrs.get(key, None)
        mock_ds = self.mock_dataset
        mock_interp_sza_s2.return_value = 1

        # Set up the mock for the S2RutAlgo class attributes
        mock_algo_instance = mock_s2_rut_algo.return_value
        mock_algo_instance.u_diff_cos = 2.0
        mock_algo_instance.u_diff_k = 3.0
        mock_algo_instance.u_ADC = 4.0
        mock_algo_instance.u_gamma = 5.0
        mock_algo_instance.k = 6.0
        mock_u_diff_temp_rate = {"Sentinel-2A": [7.0]}

        with mock.patch(
            "s2_rut_python.interface.TIME_INIT",
            {"Sentinel-2A": datetime.datetime(2015, 6, 23, 10, 00)},
        ):
            with mock.patch(
                "s2_rut_python.interface.conf.u_diff_temp_rate", mock_u_diff_temp_rate
            ):
                # Initialize the S2RUT object and run the get_band_unc_parameters method
                rut = S2RUT()
                result = rut.get_band_unc_parameters(mock_ds, "B01")

                # Define the expected keys
                expected_keys = [
                    "a",
                    "e_sun",
                    "u_sun",
                    "tecta",
                    "quant",
                    "alpha",
                    "beta",
                    "u_diff_cos",
                    "u_diff_k",
                    "u_diff_temp",
                    "u_ADC",
                    "u_gamma",
                    "k",
                ]

                # Check if all expected keys are in the result dictionary
                print(result)
                for key in expected_keys:
                    self.assertIn(key, result)

                expected_values = {
                    "a": 1.0,
                    "e_sun": 1.0,
                    "u_sun": 0.05,
                    "tecta": 1.0,
                    "quant": 10000,
                    "alpha": 1.0,
                    "beta": 1.0,
                    "u_diff_cos": 2.0,
                    "u_diff_k": 3.0,
                    "u_diff_temp": 7 * -0.06023271731690623,
                    "u_ADC": 4.0,
                    "u_gamma": 5.0,
                    "k": 6.0,
                }

                for key in expected_values:
                    self.assertEqual(result[key], expected_values[key])

    def test_return_unc_components(self):
        s2rut = S2RUT()
        unc_comp = s2rut.return_unc_components()

        # Check if the returned dictionary is structured as expected
        self.assertIn("random", unc_comp)
        self.assertIn("structured", unc_comp)
        self.assertIn("systematic", unc_comp)

    @mock.patch("s2_rut_python.interface.S2RUT.get_band_unc_parameters")
    @mock.patch("s2_rut_python.interface.MyS2RUTAlgo.unc_calculation")
    def test_run_ComponentsFalse(
        self, mock_unc_calculation, mock_get_band_unc_parameters
    ):
        # Mock the return value of the unc_calculation method
        mock_unc_calculation.return_value = np.random.rand(5, 5)
        mock_get_band_unc_parameters.return_value = {
            "a": 1.0,
            "e_sun": 1.0,
            "u_sun": 1.0,
            "tecta": np.ones(
                (10, 10)
            ),  # cannot get tecta to return as a 2d array (or any value) that isn't a class instance
            "quant": 1.0,
            "alpha": 1.0,
            "beta": 1.0,
            "u_diff_cos": 1.0,
            "u_diff_k": 1.0,
            "u_diff_temp": 7.532935146264837,
            "u_ADC": 1.0,
            "u_gamma": 1.0,
            "k": 1.0,
        }
        band_names = ["B01"]
        components = False

        result = self.s2_rut.run(
            self.mock_dataset, band_names, components
        )  # mock_dataset keeps giving attribute errors, as the xarray read into the code is very complex and has a lot of attributes, coords, values needed to run the method

        self.assertIn("B01", result)
        self.assertIn("total", result["B01"])
        self.assertEqual(result["B01"]["total"].shape, (5, 5))
        mock_get_band_unc_parameters.assert_called_once_with(self.mock_dataset, "B01")
        mock_unc_calculation.assert_called_once()

    @mock.patch("s2_rut_python.interface.S2RUT.get_band_unc_parameters")
    @mock.patch("s2_rut_python.interface.MyS2RUTAlgo.unc_calculation")
    def test_run_ComponentsTrue(
        self, mock_unc_calculation, mock_get_band_unc_parameters
    ):
        # Mock the return value of the unc_calculation method
        mock_unc_calculation.return_value = np.random.rand(5, 5)
        mock_get_band_unc_parameters.return_value = {
            "a": 1.0,
            "e_sun": 1.0,
            "u_sun": 1.0,
            "tecta": np.ones((10, 10)),
            "quant": 1.0,
            "alpha": 1.0,
            "beta": 1.0,
            "u_diff_cos": 1.0,
            "u_diff_k": 1.0,
            "u_diff_temp": 7.532935146264837,
            "u_ADC": 1.0,
            "u_gamma": 1.0,
            "k": 1.0,
        }

        band_names = ["B01"]
        components = True
        result = self.s2_rut.run(self.mock_dataset, band_names, components)

        self.assertIn("B01", result)
        self.assertIn("structured", result["B01"])
        self.assertIn("systematic", result["B01"])
        self.assertIn("random", result["B01"])
        self.assertEqual(result["B01"]["structured"].shape, (5, 5))
        self.assertEqual(result["B01"]["systematic"].shape, (5, 5))
        self.assertEqual(result["B01"]["random"].shape, (5, 5))

        np.testing.assert_equal(mock_unc_calculation.call_args[0][0], np.ones((5, 5)))
        self.assertEqual(mock_unc_calculation.call_args[0][1], 0)
        self.assertEqual(mock_unc_calculation.call_args[0][2], "Sentinel-2A")
        self.assertEqual(mock_unc_calculation.call_count, 3)

    # @mock.patch('s2_rut_python.interface.MyS2RUTAlgo.__init__', return_value=None)  # Prevent actual initialization
    # @mock.patch('s2_rut_python.interface.MyS2RUTAlgo.unc_calculation')
    # def test_run_parameter_overwrite(self, mock_unc_calculation, mock_init):
    #     # Mock the return value of the unc_calculation method
    #     mock_unc_calculation.return_value = np.random.rand(5, 5)
    #     band_names = ["B01"]
    #     components = False
    #     result = self.s2_rut.run(self.mock_dataset, band_names, components)
    #     # Verify that MyS2RUTAlgo was initialized with the correct parameters
    #     self.assertTrue(mock_init.called)
    #     init_args = mock_init.call_args[0][0]  # First argument to __init__
    #     self.assertIsInstance(init_args, dict)  # Ensure parameters are passed as a dictionary
    #     self.assertIn('a', init_args)
    #     self.assertIn('e_sun', init_args)
    #     self.assertIn('u_diff_k', init_args)
    #     self.assertIn("B01", result)
    #     self.assertIn(False, result["B01"])
    #     self.assertEqual(result["B01"].shape, (5, 5))

    @mock.patch("s2_rut_python.interface.S2RUT.get_band_unc_parameters")
    @mock.patch("s2_rut_python.interface.S2RUT.return_unc_components")
    def test_run_ComponentsInvalid(
        self, mock_return_unc_components, mock_get_band_unc_parameters
    ):
        band_names = ["B01"]
        components = "invalid"
        with self.assertRaises(ValueError):
            self.s2_rut.run(self.mock_dataset, band_names, components)
        mock_return_unc_components.assert_called_once()
        mock_get_band_unc_parameters.assert_called_once()


if __name__ == "__main__":
    unittest.main()
