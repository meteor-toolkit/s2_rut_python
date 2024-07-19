"""
Created on Thurs Oct 06 11:51:21 2022

From Fri 14 Oct 2022 sat2xr.py used to read in satellite data and metadata

@author: mattea.goalen
"""
import os
import glob
import importlib
import datetime
import rioxarray as rxr
import xarray as xr

snap_rut = importlib.import_module("snap-rut")
s2_rut_algo = __import__("snap-rut.src.main.python.s2_rut_algo")
rad_conf = __import__("snap-rut.src.main.python.s2_l1_rad_conf")
# from snap_rut.src.main.python import s2_rut_algo

try:
    import xml.etree.cElementTree as ET  # C implementation is much faster and consumes significantly less memory
except ImportError:
    import xml.etree.ElementTree as ET

S2_MSI_TYPE_STRING = "S2MSI1C"
res_10m = ['B02', 'B03', 'B04', 'B08']
res_20m = ['B05', 'B06', 'B07', 'B8A', 'B11', 'B12']
res_60m = ['B01', 'B09', 'B10']


class S2Rut:
    def __init__(self):
        self.filepath = None
        self.datastrip_meta = None
        self.product_meta = None
        self.spacecraft = None
        # possible values are "Sentinel-2A" and "Sentinel-2B". Used as a dictionary key
        self.s2_rut_info = None
        self.source_product = None
        self.mask_group = None
        self.rut_algo = s2_rut_algo.S2RutAlgo()
        self.unc_band = None
        self.toa_band = None
        # S2A launch date 23-june-2015 and S2A launch date 7-march-2017, time is indifferent.
        self.time_init = {'Sentinel-2A': datetime.datetime(2015, 6, 23, 10, 00),
                          'Sentinel-2B': datetime.datetime(2017, 3, 7, 10, 00)}
        self.sourceBandMap = None
        self.targetBandList = []
        self.inforoot = None
        self.rut_product_meta = None
        self.source_sza = None

    def process(self, filepath, bands=None, k=None, uncs=None):
        """

        :param filepath:
        :param bands:
        :param k:
        :param uncs:
        :return:
        """

        self.filepath = filepath  # make filepath variable accessible to all

        # if no bands selected default to processing all bands
        if bands is None:
            bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
            print('All bands selected')

        # if coverage factor (k) not defined set as 1
        if k is None:
            k = 1

        # if no uncertainty contributions selected default to all
        if uncs is None:
            uncs = [True, True, True, True, True, True, True, True, True, True, True, True]

        # contributors = ["INSTRUMENT_NOISE", "OOF_STRAYLIGHT-SYSTEMATIC", "OOF_STRAYLIGHT-RANDOM", "CROSSTALK",
        #                 "ADC_QUANTISATION", "DS_STABILITY", "GAMMA_KNOWLEDGE", "DIFFUSER-ABSOLUTE_KNOWLEDGE",
        #                 "DIFFUSER-TEMPORAL_KNOWLEDGE", "DIFFUSER-COSINE_EFFECT", "DIFFUSER-STRAYLIGHT_RESIDUAL",
        #                 "L1C_IMAGE_QUANTISATION"]

        self._get_metadata()

        # get filepaths to satellite data
        img_data = glob.glob(filepath + '\\*\\*\\IMG_DATA\\*B*.jp2')
        qi_detfoo = glob.glob(filepath + '\\*\\*\\QI_DATA\\*DETFOO_B*.jp2')
        qi_quality = glob.glob(filepath + '\\*\\*\\QI_DATA\\*QUALIT_B*.jp2')
        b = [img_data[i][-7:-4] for i in range(len(img_data))]

        # select only data from bands requested
        img = [img_data[i] for i in range(len(img_data)) if b[i] in bands]
        qid = [qi_detfoo[i] for i in range(len(qi_detfoo)) if b[i] in bands]
        qiq = [qi_quality[i] for i in range(len(qi_quality)) if b[i] in bands]

        # open files and set values in xarray Dataset
        self.s2_ds = xr.Dataset()
        for (bnd, im, qd, qq) in zip(bands, img, qid, qiq):
            self.s2_ds[bnd] = rxr.open_rasterio(im)[0]
            self.s2_ds[bnd + "_qiq"] = rxr.open_rasterio(qq)
            # self.s2_ds[bnd] = [rxr.open_rasterio(im), rxr.open_rasterio(qd), rxr.open_rasterio(qq)]
            # self.s2_ds[bnd][2] = self.s2_ds[bnd][2].assign_coords({'band': ["ancillary_lost", "ancillary_degraded",
            #                                                                 "msi_lost", "msi_degraded", "defective",
            #                                                                 "nodata", "partially_corrected_crosstalk",
            #                                                                 "saturated_l1a"]})
            # self.s2_ds[bnd][2] = self.s2_ds[bnd][2].rename({'band': 'mask'})

        self.toa_band_names = bands

        self.rut_algo.k = k
        self.rut_algo.unc_select = uncs
        self.source_sza = self.get_tecta()  # todo function defined below to resample sza

        # todo set attributes to add to obsarray metadata at the end
        # Source Product
        # Coverage Factor
        # Version
        # Uncertainty Contributors

        # compute uncertainty calculations
        for b in self.toa_band_names:
            self._run_s2_rut(b)

    def _get_metadata(self):
        """

        :return:
        """
        self.datastrip_meta = ET.parse(glob.glob(self.filepath + '\\*\\*\\MTD_DS.xml')[0]).getroot()
        self.product_meta = ET.parse(glob.glob(self.filepath + '\\*L1C.xml')[0]).getroot()
        self.spacecraft = self.product_meta.findall(".//SPACECRAFT_NAME")[0].text
        product_type = self.product_meta.findall(".//PRODUCT_TYPE")[0].text

        if product_type is not S2_MSI_TYPE_STRING:
            raise RuntimeError('Source product must be of type "' + S2_MSI_TYPE_STRING + '"')
        # verify as valid Sentinel-2 product first
        else:
            pass

        self.rut_algo.u_sun = self.product_meta.findall(".//U")[0].text
        self.rut_algo.quant = self.product_meta.findall(".//QUANTIFICATION_VALUE")[0].text

        # get cloud masks
        cloud_file = glob.glob(self.filepath + '\\*\\*\\QI_DATA\\*CLASSI_B*.jp2')
        self.cloud_masks = rxr.open_rasterio(cloud_file[0])
        self.cloud_masks = self.cloud_masks.assign_coords(
            {'band': ["opaque_clouds", "cirrus_clouds", "snow_and_ice_areas"]})
        self.cloud_masks = self.cloud_masks.rename({'band': 'mask'})
        # "opaque_clouds", "cirrus_clouds", "snow_and_ice_areas"

        b_10 = [res_10m[i] for i in range(len(res_10m)) if res_10m[i] in self.toa_band_names]
        b_20 = [res_20m[i] for i in range(len(res_20m)) if res_20m[i] in self.toa_band_names]
        b_60 = [res_60m[i] for i in range(len(res_60m)) if res_60m[i] in self.toa_band_names]

        if 'B02' or 'B03' or 'B04' or 'B08' in self.toa_band_names:
            self.cloud_masks_10m = self.cloud_masks.interp_like(self.s2_ds[b_10[0]])
        if 'B05' or 'B06' or 'B07' or 'B8A' or 'B11' or 'B12' in self.toa_band_names:
            self.cloud_masks_20m = self.cloud_masks.interp_like(self.s2_ds[b_20[0]])
        if 'B01' or 'B09' or 'B10' in self.toa_band_names:
            self.cloud_masks_60m = self.cloud_masks

    def _run_s2_rut(self, band):
        sampling = self.product_meta.findall(".//QUANTIFICATION_VALUE")[0].text

    def get_tecta(self):
        pass
