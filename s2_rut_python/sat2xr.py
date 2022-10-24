"""
Created Fri Oct 14 2022

@author : mg13
"""
import numpy as np
import datetime as dt
import glob
import importlib

import xarray as xr
import rioxarray as rxr

from pyproj import Transformer

snap_rut = importlib.import_module("snap-rut")
# s2_rut_algo = __import__("snap-rut.src.main.python.s2_rut_algo")
# rad_conf = __import__("snap-rut.src.main.python.s2_l1_rad_conf")

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


class Sat2xr:
    """
    Class to read in satellite image data and metadata into xarrays
    ...
    Attributes
    ----------
    filepath : str
        filepath to satellite data folder
    bands: list
        desired bands to be read in (defaults to all)

    Methods
    -------
    process(img_data, lon, lat):
        Reads data into xarray
    convert_xy(im_input, input_format, output_format="EPSG:4326"):
        Converts coordinates from one CRS to another
    clipbox_coords(lon, lat):
        Clips satellite data to limits defined by lon and lat
    """

    def __init__(self, filepath: str, bands: list = None):
        """
        Constructs all initial attributes for the satellite data object

        :param filepath: filepath to satellite data folder
        :param bands: desired bands to be read in (defaults to all)
        """
        self.filepath = filepath
        self.bands = bands
        self.lv_dict = {}
        self.gv_dict = {}
        self.msk_dict = {}
        self._res = []
        self.sat_ds = xr.Dataset()

    def process(self, img_data: list, lon: tuple or list = None, lat: tuple or list = None):
        """
        Reads in selected satellite data, populates sat_ds xarray and assigns band specific metadata attributes.

        :param img_data: list of filepaths to band image files
        :param lon: min and max longitude limits (default is None)
        :param lat: min and max latitude limits (default is None)
        :return: None
        """
        self.sat_ds = self.sat_ds.assign_attrs(self.gv_dict)
        for bnd in self.bands:
            self.sat_ds[bnd] = None
        res = [*set(self._res)]
        for i in res:
            for (bnd, im, r) in zip(self.bands, img_data, self._res):
                if r == i:
                    if lat is not None:
                        pl_hol = rxr.open_rasterio(im)[0].rio.clip_box(lon[0], lat[0], lon[1], lat[1], crs="EPSG:4326")
                    else:
                        pl_hol = rxr.open_rasterio(im)[0]
                    lon_new, lat_new = self.convert_xy((pl_hol.x, pl_hol.y), pl_hol.rio.crs, "EPSG:4326")
                    self.lv_dict[bnd].update({'LONGITUDE': lon_new, 'LATITUDE': lat_new})
                    pl_hol = pl_hol.rename({'x': 'x_' + r, 'y': 'y_' + r})
                    pl_hol = pl_hol.assign_attrs(self.lv_dict[bnd])
                    self.sat_ds[bnd] = pl_hol
                else:
                    pass

    @staticmethod
    def convert_xy(im_input: list or tuple or np.ndarray, input_format: str, output_format: str = "EPSG:4326"):
        """
        Convert input coordinates from their coordinate reference system to the World Geodetic System 1984 (WGS 84) or
        other specified coordinate reference system.

        If im_input is a list/1D- array of x and y coordinates instead of a two dimensional coordinate grid, list/array
        first used to create a two dimensional coordinate grid

        :param im_input: list/1D-array of coordinate values, or 2D-coordinate grid
        :param input_format: coordinate reference system of input coordinates
        :param output_format: desired coordinate reference system (defaults to 'EPSG:4326')
        :return: converted coordinates
        """
        if im_input[0].shape == im_input[1].shape and im_input[0].ndim == 2:
            x_mesh = im_input[0]
            y_mesh = im_input[1]
        else:
            x_mesh, y_mesh = np.meshgrid(im_input[0], im_input[1])
        transformer = Transformer.from_crs(input_format, output_format, always_xy=True)
        out_x, out_y = transformer.transform(x_mesh, y_mesh)
        return out_x, out_y

    def clipbox_coords(self, lon: tuple or list, lat: tuple or list):  # checked by band resolution not enough memory
        """
        Clips the satellite data to desired longitude and latitude limits post file reading and saves as an object
        variable.

        :param lon: min and max longitude limits
        :param lat: min and max latitude limits
        :return: None
        """
        clip_ds = []
        for bnd in self.sat_ds:
            res = self.sat_ds[bnd].RESOLUTION
            clip_da = self.sat_ds[bnd].rename({'x_' + res: 'x', 'y_' + res: 'y'})
            clip_da = clip_da.rio.clip_box(lon[0], lat[0], lon[1], lat[1], crs="EPSG:4326")
            lon_new, lat_new = self.convert_xy((clip_da.x, clip_da.y), clip_da.rio.crs, "EPSG:4326")
            self.lv_dict[bnd].update({'LONGITUDE': lon_new, 'LATITUDE': lat_new})
            clip_da = clip_da.rename({'x': 'x_' + res, 'y': 'y_' + res})
            clip_da = clip_da.assign_attrs({'LONGITUDE': lon_new, 'LATITUDE': lat_new})
            clip_ds.append(clip_da.to_dataset())
        self.clip_ds = xr.merge(clip_ds)
        self.clip_ds = self.clip_ds.assign_attrs(self.gv_dict)


class S2toxr(Sat2xr):
    """
    Subclass of Sat2xr to read from Sentinel-2A and 2B
    ...
    Attributes
    ----------
    filepath : str
        filepath to satellite data folder
    bands: list
        desired bands to be read in (defaults to all)

    Methods
    -------
    read_meta():
        Reads in metadata and assigns to object variable.
    read_img(lon: tuple or list = None, lat: tuple or list = None):
        Reads satellite data into sat_ds xarray and assigns band specific attributes.
    get_masks():
        Gets quality and cloud mask variables and saves to object.
    """

    def __init__(self, filepath: str, bands: list = None):  # define desired variables
        """
        Constructs all initial attributes for the Sentinel-2 data object.

        :param filepath: filepath to satellite data folder
        :param bands: desired bands to be read in (defaults to all)
        """
        super().__init__(filepath, bands)
        self._bands_dict = {'B01': '0', 'B02': '1', 'B03': '2', 'B04': '3', 'B05': '4', 'B06': '5', 'B07': '6',
                            'B08': '7', 'B8A': '8', 'B09': '9', 'B10': '10', 'B11': '11', 'B12': '12'}
        self._SOLAR_IRRADIANCE_TAG = 'SOLAR_IRRADIANCE'
        self._RESOLUTION_TAG = 'RESOLUTION'
        self._ALPHA_TAG = 'ALPHA'
        self._BETA_TAG = 'BETA'
        self._PHYSICAL_GAINS_TAG = 'PHYSICAL_GAINS'
        self._SPACECRAFT_NAME_TAG = 'SPACECRAFT_NAME'
        self._DATATAKE_SENSING_START_TAG = 'DATATAKE_SENSING_START'
        self._QUANTIFICATION_VALUE_TAG = 'QUANTIFICATION_VALUE'
        self._U_TAG = 'U'
        self._TIME_DEGRADATION_TAG = 'TIME_DEGRADATION'
        self._time_init = {'Sentinel-2A': dt.datetime(2015, 6, 23, 10, 00),
                           'Sentinel-2B': dt.datetime(2017, 3, 7, 10, 00)}
        self._u_diff_temp_rate = {
            'Sentinel-2A': [0.15, 0.09, 0.04, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'Sentinel-2B': [0.15, 0.09, 0.04, 0.02, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}

    def read_meta(self, sza=False):
        """
        Reads in the selected metadata and sets attributes to object global and band specific dictionary.

        If no bands are selected prints 'All bands selected' and reads all band metadata.

        :param sza: whether to read in the solar zenith angles or not
        :return:
        """
        try:
            datastrip_meta = ET.parse(glob.glob(self.filepath + '\\*\\*\\MTD_DS.xml')[0]).getroot()
            granule_meta = ET.parse(glob.glob(self.filepath + '\\*\\*\\MTD_TL.xml')[0]).getroot()
        except IndexError:
            raise RuntimeError('Source product metadata not found in "' + self.filepath + '"')

        if self.bands is None:
            print('All bands selected')
            self.bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
        else:
            pass

        # set global variables
        self._spacecraft = datastrip_meta.findall(".//" + self._SPACECRAFT_NAME_TAG)[0].text
        self._data_ss = dt.datetime.strptime(datastrip_meta.findall(".//" + self._DATATAKE_SENSING_START_TAG)[0].text,
                                             '%Y-%m-%dT%H:%M:%S.%fZ')
        quant = datastrip_meta.findall(".//" + self._QUANTIFICATION_VALUE_TAG)[0].text
        u = datastrip_meta.findall(".//U")[0].text

        self.gv_dict = {self._SPACECRAFT_NAME_TAG: self._spacecraft, self._QUANTIFICATION_VALUE_TAG: quant,
                        self._U_TAG: u}

        # set local variables
        no_d = (self._data_ss - self._time_init[self._spacecraft]).days / 365.25

        for bnd in self.bands:
            self.lv_dict[bnd] = {
                self._ALPHA_TAG: datastrip_meta.findall(".//*[@bandId='" + self._bands_dict[bnd] + "']/*/ALPHA")[
                    0].text,
                self._RESOLUTION_TAG:
                    datastrip_meta.findall(".//*[@bandId='" + self._bands_dict[bnd] + "']/RESOLUTION")[0].text,
                self._BETA_TAG: datastrip_meta.findall(".//*[@bandId='" + self._bands_dict[bnd] + "']/*/BETA")[
                    0].text,
                self._PHYSICAL_GAINS_TAG:
                    datastrip_meta.findall(".//*[@bandId='" + self._bands_dict[bnd] + "']/PHYSICAL_GAINS")[0].text,
                self._SOLAR_IRRADIANCE_TAG:
                    datastrip_meta.findall(".//SOLAR_IRRADIANCE[@bandId='" + self._bands_dict[bnd] + "']")[0].text,
                self._TIME_DEGRADATION_TAG: no_d * self._u_diff_temp_rate[self._spacecraft][int(self._bands_dict[bnd])]}
            self._res.append(
                datastrip_meta.findall(".//*[@bandId='" + self._bands_dict[bnd] + "']/RESOLUTION")[0].text)

        if sza == True:
            self.read_sza(granule_meta)
        else:
            pass

    def read_img(self, lon: tuple or list = None, lat: tuple or list = None):
        """
        Finds filepaths to desired image files, reads in satellite data into sat_ds xarray and assigns band specific
        attributes.

        :param lon: min and max longitude limits
        :param lat: min and max latitude limits
        :return: None
        """
        img_data = [glob.glob(self.filepath + '\\*\\*\\IMG_DATA\\*' + i + '.jp2')[0] for i in self.bands]
        # populate Dataset and assign attributes
        super().process(img_data, lon, lat)

        # return self.sat_ds

    def get_masks(self):
        """
        Reads, interpolates and saves quality and cloud mask variables to object in band specific dictionary and as
        band specific attributes.

        :return: None
        """
        if self._data_ss > dt.datetime(2021, 10, 26):
            msk_quality = [glob.glob(self.filepath + '\\*\\*\\QI_DATA\\*QUALIT_' + i + '.jp2')[0] for i in self.bands]
            # detector footprints for each band
            msk_dfoot = [glob.glob(self.filepath + '\\*\\*\\QI_DATA\\*DETFOO_' + i + '.jp2')[0] for i in self.bands]
            msk_clouds = glob.glob(self.filepath + '\\*\\*\\QI_DATA\\*CLASSI*.jp2')[0]

            xds_clouds = rxr.open_rasterio(msk_clouds)

            clouds_dict = {}

            for i, bnd in enumerate(self.bands):
                self.msk_dict[bnd] = {"detector_footprints": rxr.open_rasterio(msk_dfoot[i])[0],
                                      "ancillary_lost": rxr.open_rasterio(msk_quality[i])[0],
                                      "ancillary_degraded": rxr.open_rasterio(msk_quality[i])[1],
                                      "msi_lost": rxr.open_rasterio(msk_quality[i])[2],
                                      "msi_degraded": rxr.open_rasterio(msk_quality[i])[3],
                                      "defective": rxr.open_rasterio(msk_quality[i])[4],
                                      "nodata": rxr.open_rasterio(msk_quality[i])[5],
                                      "partially_corrected_crosstalk": rxr.open_rasterio(msk_quality[i])[6],
                                      "saturated_l1a": rxr.open_rasterio(msk_quality[i])[7]}
                try:
                    self.msk_dict[bnd].update(
                        {"opaque_clouds": clouds_dict[self._res[i]][0], "cirrus_clouds": clouds_dict[self._res[i]][1],
                         "snow_and_ice_areas": clouds_dict[self._res[i]][2]})
                except KeyError:
                    clouds_dict[self._res[i]] = xds_clouds.rio.reproject_match(
                        self.sat_ds[bnd].rename({'x_' + self._res[i]: 'x', 'y_' + self._res[i]: 'y'}))
                    self.msk_dict[bnd].update(
                        {"opaque_clouds": clouds_dict[self._res[i]][0], "cirrus_clouds": clouds_dict[self._res[i]][1],
                         "snow_and_ice_areas": clouds_dict[self._res[i]][2]})

    def read_sza(self, granule_meta):
        """
        Reads sza Viewing Incidence Grid for desired bands

        :param granule_meta:
        :return:
        """
        for bnd in self.bands:
            test_sza = []
            for dec in range(5):
                test = [i.text.split() for i in granule_meta.findall(
                    ".//Viewing_Incidence_Angles_Grids[@bandId='" + self._bands_dict[bnd] + "'][@detectorId='" + str(
                        dec + 1) + "']/Zenith/Values_List/VALUES")]
                for j in range(len(test[0])):
                    test[j] = [float(k) if float(k) > 0 else 0 for k in test[j]]
                test_sza.append(test)
            test_sza_np = np.sum(np.array(test_sza), axis=0)

            msk_sza_test = []
            for dec in range(5):
                msk_test = []
                for i in range(23):
                    msk_test.append([1 if float(k) > 0 else 0 for k in test_sza[dec][i]])
                msk_sza_test.append(msk_test)
            msk_sza = np.sum(np.array(msk_sza_test), axis=0)

            sza_final = np.nan_to_num(np.divide(test_sza_np, msk_sza), nan=0)
            self.lv_dict[bnd].update({'SOLAR_ZENITH_ANGLES': sza_final})
