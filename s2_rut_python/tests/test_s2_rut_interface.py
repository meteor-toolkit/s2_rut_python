import datetime
import os, sys
import unittest
import unittest.mock as mock
from unittest.mock import MagicMock
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
                #     'tecta': np.ones((10, 10)),
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

    @mock.patch('s2_rut_python.s2_rut_interface.MyS2RUTAlgo')  # Mock MyS2RUTAlgo class
    @mock.patch('s2_rut_python.s2_rut_interface.S2RUT.return_unc_correlations')  # Mock return_unc_correlations method
    @mock.patch('s2_rut_python.s2_rut_interface.S2RUT.get_band_unc_parameters')  # Mock get_band_unc_parameters method
    def test_run_rut(self, mock_get_band_unc_parameters, mock_return_unc_correlations, mock_my_s2_rut_algo):


        # Setup the mock return values
        mock_return_unc_correlations.return_value = {
            'systematic': {'1': True, '2': False},
            'random': {'1': True, '2': False},
            'structured': {'1': True, '2': False}
        }

        mock_get_band_unc_parameters.return_value = {
            'a': 1.0,
            'e_sun': 1.0,
            'u_sun': 1.0,
            'tecta': np.ones((10980, 10980)),
            'quant': 1.0,
            'alpha': 1.0,
            'beta': 1.0,
            'u_diff_cos': 1.0,
            'u_diff_k': 1.0,
            'u_diff_temp': 7.532935146264837,
            'u_ADC': 1.0,
            'u_gamma': 1.0,
            'k': 1.0
        }

        mock_rut = mock_my_s2_rut_algo.return_value
        mock_rut.unc_select = [True, True, True, True, True, True, True, True, True, True, True, True, True]
        mock_rut.unc_calculation.return_value = np.ones((10, 10))

        # Create a mock xarray.Dataset
        mock_dataset = xr.Dataset({
            'meas_vars': 'B01',
            'B01': (('x', 'y'), np.random.rand(10, 10)),
            'platform': 'Sentinel-2A'
        })

        # Initialize the S2RUT object
        s2rut = S2RUT()

        mock_rut.unc_calculation(mock_dataset, 0, 'Sentinel-2A')

        # Test with 'components' option
        result_comp = s2rut.run_rut(mock_dataset, band_names=['B01'], unc_info='components')
        self.assertIn('B01', result_comp)
        self.assertIn('structured', result_comp['B01'])
        self.assertIn('systematic', result_comp['B01'])
        self.assertIn('random', result_comp['B01'])
        self.assertIsInstance(result_comp['B01']['systematic'], xr.DataArray)

        # Test with 'total' option
        result_total = s2rut.run_rut(mock_dataset, band_names=['B01'], unc_info='total')
        self.assertIn('B01', result_total)
        self.assertIn('total', result_total['B01'])
        self.assertIsInstance(result_total['B01']['total'], xr.DataArray)

        # Test with invalid 'unc_info'
        with self.assertRaises(ValueError):
            s2rut.run_rut(mock_dataset, band_names=['B01'], unc_info='invalid_option')


if __name__ == "__main__":
    unittest.main()
