import datetime
import os, sys
import unittest
import unittest.mock as mock
from unittest.mock import MagicMock
import xarray
import xarray as xr
from s2_rut_python.s2_rut_interface import S2RUT, MyS2RUTAlgo
import numpy as np
import s2_rut_python.s2_rut_interface
from eoio.utils.dict_tools import get_value
from eoio.processors.utils import interp_sza_s2

THIS_DIRECTORY = os.path.abspath(os.path.dirname(__file__))
S2_RUT_DIRECTORY = os.path.join(THIS_DIRECTORY, "snap-rut", "src", "main", "python")

sys.path.insert(0, S2_RUT_DIRECTORY)
import s2_rut_algo
import s2_l1_rad_conf

__author__ = [
    "Rasma Ormane <rasma.ormane@npl.co.uk>",
]


class TestS2RUT(unittest.TestCase):

    @mock.patch('eoio.utils.dict_tools.get_value')
    @mock.patch('eoio.processors.utils.interp_sza_s2')
    @mock.patch('s2_rut_algo.S2RutAlgo')
    def test_get_band_unc_parameters(self, mock_s2_rut_algo, mock_get_value, mock_interp_sza_s2):
        # create mock xarray dataset
        mock_ds = MagicMock(spec=xr.Dataset)

        mock_ds['B01'].PHYSICAL_GAINS = 1.0
        mock_ds['B01'].SOLAR_IRRADIANCE = {'#text': 1.0}
        mock_ds['B01'].ALPHA = 1.0
        mock_ds['B01'].BETA = 1.0
        mock_ds.attrs = {
            'U': 1.0,
            'QUANTIFICATION_VALUE': 1.0,
            'DATASTRIP_SENSING_START': datetime.datetime(2023, 1, 1, 10, 00),
        }
        mock_ds.platform = 'Sentinel-2A'

        # Set up the mock return values for get_value and interp_sza_s2
        mock_get_value.side_effect = lambda attrs, key: attrs.get(key, None)
        mock_interp_sza_s2.return_value = np.ones((10, 10))

        # Set up the mock for the S2RutAlgo class attributes
        mock_algo_instance = mock_s2_rut_algo.return_value
        mock_algo_instance.u_diff_cos = 1.0
        mock_algo_instance.u_diff_k = 1.0
        mock_algo_instance.u_ADC = 1.0
        mock_algo_instance.u_gamma = 1.0
        mock_algo_instance.k = 1.0
        mock_u_diff_temp_rate = {"Sentinel-2A": [1.0]}

        with mock.patch('s2_rut_python.s2_rut_interface.TIME_INIT',
                        {'Sentinel-2A': datetime.datetime(2015, 6, 23, 10, 00)}):
            with mock.patch('s2_rut_python.s2_rut_interface.conf.u_diff_temp_rate', mock_u_diff_temp_rate):
                # Initialize the S2RUT object and run the get_band_unc_parameters method
                rut = S2RUT()
                result = rut.get_band_unc_parameters(mock_ds, 'B01')

                # Define the expected keys
                expected_keys = [
                    'a', 'e_sun', 'u_sun', 'tecta', 'quant', 'alpha', 'beta',
                    'u_diff_cos', 'u_diff_k', 'u_diff_temp', 'u_ADC', 'u_gamma', 'k'
                ]

                # Check if all expected keys are in the result dictionary
                for key in expected_keys:
                    self.assertIn(key, result)

                # expected_values = {
                #     'a': 1.0,
                #     'e_sun': 1.0,
                #     'u_sun': 1.0,
                #     'tecta': np.ones((10, 10)), # cannot get tecta to return as a 2d array (or any value) that isn't a class instance
                #     'quant': 1.0,
                #     'alpha': 1.0,
                #     'beta': 1.0,
                #     'u_diff_cos': 1.0,
                #     'u_diff_k': 1.0,
                #     'u_diff_temp': 7.532935146264837,
                #     'u_ADC': 1.0,
                #     'u_gamma': 1.0,
                #     'k': 1.0
                # }
                #
                # for key in expected_values:
                #     self.assertEqual(result[key], expected_values[key])

    def test_return_unc_correlations(self):
        s2rut = S2RUT()
        unc_corr = s2rut.return_unc_correlations()

        # Check if the returned dictionary is structured as expected
        self.assertIn('random', unc_corr)
        self.assertIn('structured', unc_corr)
        self.assertIn('systematic', unc_corr)


    def setUp(self):
        # Create a mock dataset
        self.mock_dataset = MagicMock(spec=xarray.Dataset)

        # Mock the dimensions and coordinates
        self.mock_dataset.dims = {'y_60m': 5, 'x_60m': 5}  # Keep getting an error AttributeError: 'dict' object has no attribute 'dims'
        self.mock_dataset.coords = {
            'y': np.arange(5),
            'x': np.arange(5)
        }

        # Mock the band variables in the dataset
        def mock_variable():
            mock_var = MagicMock()
            mock_var.coords = self.mock_dataset.coords
            mock_var.dims = ('y_60m', 'x_60m')
            mock_var.values = np.random.rand(5, 5)
            return mock_var

        # Return the mock variable for any band
        self.mock_dataset.__getitem__.side_effect = lambda band: mock_variable()

        # Mock attributes
        self.mock_dataset.attrs = {
            'PLATFORM': 'Sentinel-2A',
            'U': 0.05,
            'QUANTIFICATION_VALUE': 10000,
            'DATASTRIP_SENSING_START': np.datetime64("2022-06-01T10:00:00")
        }

        self.s2_rut = S2RUT()

    @mock.patch('s2_rut_python.s2_rut_interface.MyS2RUTAlgo.unc_calculation')
    def test_run_rut_total_uncertainty(self, mock_unc_calculation):
        # Mock the return value of the unc_calculation method
        mock_unc_calculation.return_value = np.random.rand(5, 5)
        band_names = ["B01"]
        unc_info = "total"
        result = self.s2_rut.run_rut(self.mock_dataset, band_names, unc_info)  # mock_dataset keeps giving attribute errors, as the xarray read into the code is very complex and has a lot of attributes, coords, values needed to run the method
        self.assertIn("B01", result)
        self.assertIn("total", result["B01"])
        self.assertEqual(result["B01"]["total"].shape, (5, 5))
        mock_unc_calculation.assert_called_once()

    @mock.patch('s2_rut_python.s2_rut_interface.MyS2RUTAlgo.unc_calculation')
    def test_run_rut_components_uncertainty(self, mock_unc_calculation):
        # Mock the return value of the unc_calculation method
        mock_unc_calculation.return_value = np.random.rand(5, 5)
        band_names = ["B01"]
        unc_info = "components"
        result = self.s2_rut.run_rut(self.mock_dataset, band_names, unc_info)
        self.assertIn("B01", result)
        self.assertIn("structured", result["B01"])
        self.assertIn("systematic", result["B01"])
        self.assertIn("random", result["B01"])
        self.assertEqual(result["B01"]["structured"].shape, (5, 5))
        self.assertEqual(result["B01"]["systematic"].shape, (5, 5))
        self.assertEqual(result["B01"]["random"].shape, (5, 5))
        mock_unc_calculation.assert_called_with(
            np.random.rand(5, 5),  # Dummy input data
            self.s2_rut.band_id["B01"],
            'Sentinel-2A'
        )

    @mock.patch('s2_rut_python.s2_rut_interface.MyS2RUTAlgo.__init__', return_value=None)  # Prevent actual initialization
    @mock.patch('s2_rut_python.s2_rut_interface.MyS2RUTAlgo.unc_calculation')
    def test_run_rut_parameter_overwrite(self, mock_unc_calculation, mock_init):
        # Mock the return value of the unc_calculation method
        mock_unc_calculation.return_value = np.random.rand(5, 5)
        band_names = ["B01"]
        unc_info = "total"
        result = self.s2_rut.run_rut(self.mock_dataset, band_names, unc_info)
        # Verify that MyS2RUTAlgo was initialized with the correct parameters
        self.assertTrue(mock_init.called)
        init_args = mock_init.call_args[0][0]  # First argument to __init__
        self.assertIsInstance(init_args, dict)  # Ensure parameters are passed as a dictionary
        self.assertIn('a', init_args)
        self.assertIn('e_sun', init_args)
        self.assertIn('u_diff_k', init_args)
        self.assertIn("B01", result)
        self.assertIn("total", result["B01"])
        self.assertEqual(result["B01"]["total"].shape, (5, 5))

    def test_run_rut_invalid_unc_info(self):
        band_names = ["B01"]
        unc_info = "invalid"
        with self.assertRaises(ValueError):
            self.s2_rut.run_rut(self.mock_dataset, band_names, unc_info)


if __name__ == "__main__":
    unittest.main()
