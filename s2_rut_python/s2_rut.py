# -*- coding: utf-8 -*-
"""
EOpy adaptation of the Sentinel-2 Radiometric Uncertainty tool

Based on the python code found in the git repository:
https://github.com/senbox-org/snap-rut/tree/master/src/main/python

- I've not changed any of the calculations
- I've only made it eopy compatible
- I've got rid of the Java stuff to do with putting the product into snappy

Created 04/03/2020 by N. Origo (original code by Javier Gorrono)
"""
import snappy
from snappy import HashMap as hash
import s2_rut_python.s2_rut_algo as s2_rut_algo
import numpy as np
import datetime
import os
import pdb

try:
    import xml.etree.cElementTree as ET  # C implementation is much faster and consumes significantly less memory
except ImportError:
    import xml.etree.ElementTree as ET
import s2_rut_python.s2_l1_rad_conf as rad_conf

# necessary for logging
# from snappy import SystemUtils

S2_MSI_TYPE_STRING = "S2_MSI_Level-1C"
S2_BAND_NAMES = [
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B8A",
    "B9",
    "B10",
    "B11",
    "B12",
]
S2_BAND_SAMPLING = {
    "B1": 60,
    "B2": 10,
    "B3": 10,
    "B4": 10,
    "B5": 20,
    "B6": 20,
    "B7": 20,
    "B8": 10,
    "B8A": 20,
    "B9": 60,
    "B10": 60,
    "B11": 20,
    "B12": 20,
}

# If a Java type is needed which is not imported by snappy by default it can be retrieved manually.
# First import jpy and then the type to be imported
from snappy import jpy

MetadataElement = jpy.get_type("org.esa.snap.core.datamodel.MetadataElement")
MetadataAttribute = jpy.get_type("org.esa.snap.core.datamodel.MetadataAttribute")


class S2RutOp:
    def __init__(self, product, params):
        self.s2_rut_info = None
        self.mask_group = None
        self.product_meta = None
        self.datastrip_meta = None
        self.spacecraft = None  # possible values are "Sentinel-2A" and "Sentinel-2B". Used as a dictionary key
        self.rut_algo = s2_rut_algo.S2RutAlgo()
        self.unc_band = None
        self.toa_band = None

        # S2A launch date 23-june-2015 and S2A launch date 7-march-2017, time is indifferent.
        self.time_init = {
            "Sentinel-2A": datetime.datetime(2015, 6, 23, 10, 00),
            "Sentinel-2B": datetime.datetime(2017, 3, 7, 10, 00),
        }
        self.sourceBandMap = None
        self.targetBandList = []
        self.inforoot = None
        self.rut_product_meta = None
        self.source_sza = None
        self.ubands = []

        try:
            resolution = float(params["resolution"])
        except:
            raise IOError(
                "You must set the resolution required as an number (generally: 10, 20 or 60)"
            )
        try:
            bname = params["band_name"]
        except:
            raise IOError("You must select a band to process")
        """
        if "band_names" in params.keys():
            self.toa_band_names = params["band_names"]
        else:
            self.toa_band_names = S2_BAND_NAMES
        """
        # make params visible to the other functions
        self.params = params
        self.source_product = product

        # in some weird case where you might have float resolutions??
        if resolution % 1 == 0:
            resolution = int(resolution)
        self.resolution = resolution

        # for resid in np.arange(len(self.product.product)):
        #     if "opaque_clouds_"+str(resolution)+"m" in product.product[resid]["variables"]:
        #         self.source_product = product.product[resid]["product"]
        #         break

        self.initialise()
        band = self.source_product.getBand(bname)
        uband, bw, bh = self.computeBand(band)
        self.ubandrs = uband.reshape(bh, bw)
        # add band to the eopy product
        # something = eopy.addBand # Sam however you want to do this

        # add metadata to the product (not sure if you want to save all of these things...
        self.rut_algo.k  # coverage factor set
        self.rut_algo.e_sun  # extraterrestrial solar spectrum
        self.rut_algo.u_sun  # presumably uncertainty in solar spectrum used
        self.rut_algo.tecta  # SZA # probabl don't want to include this
        self.rut_algo.quant  # Image quantisation prsumably
        self.rut_algo.unc_select  # uncertainty contributors selected
        self.rut_algo.a  # presumably the absolute calibration coefficient
        self.rut_algo.alpha  # noise model parameter
        self.rut_algo.beta  # noise model parameter
        self.rut_algo.u_diff_temp  # Diffuser reflectance temporal knowledge uncertainty
        self.rut_algo.u_diff_cos  # Angular diffuser knowledge uncertainty (cosine effect)
        self.rut_algo.u_diff_k  # Residual traylight in calibration mode uncertainty
        self.rut_algo.u_diff_abs  # Diffuser reflectance absolute knowledge uncertainty
        self.rut_algo.u_ADC  # ADC quantisation uncertainty
        self.rut_algo.u_gamma  # Non-linearity and non-uniformity knowledge uncertainty
        self.rut_algo.u_noise  # Instrument noise uncertainty
        self.rut_algo.u_stray_sys  # OOF stray light uncertainty - systematic part
        self.rut_algo.u_stray_rand  # OOF stray light uncertainty - random part
        self.rut_algo.u_xtalk  # Crosstalk uncertainty
        self.rut_algo.u_DS  # Dark signal stability uncertainty
        self.rut_algo.u_ref_quant  # L1C image quantisation uncertainty

        self.args  # potential args
        params.keys()  # actual args selected

        # not sure what else you'd want to put in here

        """
        ##### compute for all resolutions #####
        for p in [range(3)[-1]]:
            self.source_product   = product.product[p]["product"]
            self.initialise()
            ubands = []

            for tb in self.toa_band_names:
                band        = self.source_product.getBand(tb)
                uband,bw,bh = self.computeBand(band)
                ubandrs     = uband.reshape(bh,bw)

                # add band to the eopy product
                #something = eopy.addBand # Sam however you want to do this
                ubands.append(ubandrs)

                # add metadata to the product
                self.rut_algo.k # coverage factor set
                self.rut_algo.u_sun #
                self.rut_algo.quant #
                self.rut_algo.unc_select #
                self.rut_algo.a #
                self.rut_algo.e_sun #
                self.rut_algo.alpha #
                self.rut_algo.beta #
                self.rut_algo.u_diff_temp #
                self.args # potential args
                params.keys() # actual args
                # not sure what else you'd want to put in here

                print p,tb
            self.ubands.append(ubands)
        """

    def initialise(self):

        if self.source_product.getProductType() != S2_MSI_TYPE_STRING:
            raise RuntimeError(
                'Source product must be of type "' + S2_MSI_TYPE_STRING + '"'
            )

        self.mask_group = (
            self.source_product.getMaskGroup()
        )  # obtain the masks from the product
        metadata_root = self.source_product.getMetadataRoot()
        self.product_meta = metadata_root.getElement("Level-1C_User_Product")
        self.datastrip_meta = metadata_root.getElement("Level-1C_DataStrip_ID")

        self.spacecraft = (
            self.datastrip_meta.getElement("General_Info")
            .getElement("Datatake_Info")
            .getAttributeString("SPACECRAFT_NAME")
        )

        self.rut_algo.u_sun = self.get_u_sun(self.product_meta)
        self.rut_algo.quant = self.get_quant(self.product_meta)

        self.source_sza = self.source_product.getBand("sun_zenith")

        self.rut_algo.k = self.get_k()
        self.rut_algo.unc_select = self.get_unc_select()

    def computeBand(self, band):
        # Logging template
        # SystemUtils.LOG.info('target band name: ' + band.getName())
        # SystemUtils.LOG.info('tile rect: ' + tile.getRectangle().toString())
        source_band = band
        gname = source_band.getName()
        toa_band_id = np.int(S2_BAND_NAMES.index(gname))

        source_sza = self.source_sza

        bw = source_band.getRasterWidth()
        bh = source_band.getRasterHeight()

        # cloudmask  = self.mask_roi('opaque_clouds_'+ str(self.resolution) +"m",
        #                            0,0,bw,bh)
        # cirrusmask = self.mask_roi('cirrus_clouds_' + str(self.resolution) +"m",
        #                            0,0,bw,bh)

        szah = source_sza.getRasterHeight()
        szaw = source_sza.getRasterWidth()

        # remove this after testing
        try:
            assert szah == bh
        except:
            raise IOError("Heigth of band and solar zenith arrays are not the same")
        try:
            assert szaw == bw
        except:
            raise IOError("Width of band and solar zenith arrays are not the same")

        sza_samples = source_sza.readPixels(0, 0, szaw, szah, np.zeros(szaw * szah))
        self.rut_algo.tecta = sza_samples

        self.rut_algo.a = self.get_a(self.datastrip_meta, toa_band_id)
        self.rut_algo.e_sun = self.get_e_sun(self.product_meta, toa_band_id)
        self.rut_algo.alpha = self.get_alpha(self.datastrip_meta, toa_band_id)
        self.rut_algo.beta = self.get_beta(self.datastrip_meta, toa_band_id)
        self.rut_algo.u_diff_temp = self.get_u_diff_temp(
            self.datastrip_meta, toa_band_id
        )

        toa_samples = source_band.readPixels(0, 0, bw, bh, np.zeros(bw * bh))

        # this is the core where the uncertainty calculation should grow
        unc = self.rut_algo.unc_calculation(
            np.array(toa_samples, dtype=np.float64), toa_band_id, self.spacecraft
        )

        return unc, bw, bh

    def dispose(self, context):
        pass

    def get_quant(self, product_meta):
        return (
            product_meta.getElement("General_info")
            .getElement("Product_Image_Characteristics")
            .getAttributeDouble("QUANTIFICATION_VALUE")
        )

    def get_u_sun(self, product_meta):
        return (
            product_meta.getElement("General_Info")
            .getElement("Product_Image_Characteristics")
            .getElement("Reflectance_Conversion")
            .getAttributeDouble("U")
        )

    def get_e_sun(self, product_meta, band_id):
        return float(
            [
                i
                for i in product_meta.getElement("General_Info")
                .getElement("Product_Image_Characteristics")
                .getElement("Reflectance_Conversion")
                .getElement("Solar_Irradiance_list")
                .getAttributes()
                if i.getName() == "SOLAR_IRRADIANCE"
            ][band_id]
            .getData()
            .getElemString()
        )

    def get_u_diff_temp(self, datastrip_meta, band_id):
        # START or STOP time has no effect. We provide a degradation based on MERIS year rates
        time_start = datetime.datetime.strptime(
            datastrip_meta.getElement("General_Info")
            .getElement("Datastrip_Time_Info")
            .getAttributeString("DATASTRIP_SENSING_START"),
            "%Y-%m-%dT%H:%M:%S.%fZ",
        )
        return (
            (time_start - self.time_init[self.spacecraft]).days
            / 365.25
            * rad_conf.u_diff_temp_rate[self.spacecraft][band_id]
        )

    def get_beta(self, datastrip_meta, band_id):
        return (
            [
                i
                for i in datastrip_meta.getElement("Quality_Indicators_Info")
                .getElement("Radiometric_Info")
                .getElement("Radiometric_Quality_list")
                .getElements()
                if i.getName() == "Radiometric_Quality"
            ][band_id]
            .getElement("Noise_Model")
            .getAttributeDouble("BETA")
        )

    def get_alpha(self, datastrip_meta, band_id):
        return (
            [
                i
                for i in datastrip_meta.getElement("Quality_Indicators_Info")
                .getElement("Radiometric_Info")
                .getElement("Radiometric_Quality_list")
                .getElements()
                if i.getName() == "Radiometric_Quality"
            ][band_id]
            .getElement("Noise_Model")
            .getAttributeDouble("ALPHA")
        )

    def get_a(self, datastrip_meta, band_id):
        return [
            i
            for i in datastrip_meta.getElement("Image_Data_Info")
            .getElement("Sensor_Configuration")
            .getElement("Acquisition_Configuration")
            .getElement("Spectral_Band_Info")
            .getElements()
            if i.getName() == "Spectral_Band_Information"
        ][band_id].getAttributeDouble("PHYSICAL_GAINS")

    def get_k(self):
        """ Assumes k=1 if nothing specified """
        if "coverage_factor" in self.params.keys():
            return self.params["coverage_factor"]
        else:
            return 1.0

    def get_unc_select(self):
        """ Parse the params """

        # get param keys
        pkeys = self.params.keys()

        # potential args
        args = [
            "ADC_quantisation",
            "Crosstalk",
            "Diffuser-absolute_knowledge",
            "Diffuser-cosine_effect",
            "Diffuser-straylight_residual",
            "Diffuser-temporal_knowledge",
            "DS_stability",
            "Gamma_knowledge",
            "Instrument_noise",
            "L1C_image_quantisation",
            "OOF_straylight-random",
            "OOF_straylight-systematic",
        ]
        self.args = args

        # all contributions flag
        if "all_contribs" in pkeys:
            if self.params["all_contribs"]:
                return np.ones(12).astype(bool)
            else:
                return np.zeros(12).astype(bool)

        else:
            # otherwise sort out the individual settings
            settings = []
            for arg in args:
                if arg in pkeys:
                    if self.params[arg]:
                        settings.append(True)
                    else:
                        settings.append(False)
                else:
                    # if setting not in there then set to False as default
                    settings.append(False)
            return np.array(settings)

    def get_masterband(self, targetBandList):
        max_width = -1
        band_index = -1
        for index, band in enumerate(targetBandList):
            width = band.getRasterWidth()
            if width > max_width:
                band_index = index
                max_width = width
        return targetBandList[band_index]

    def mask_roi(self, masktag, x, y, width, height):
        """
        The function supports the automatic read of ROI masks in function masks_extract
        :param masktag: the tag of the mask from the S2 L1C product (list of them in self.source_product.getBandNames())
        :return: ROI of raster data from the specific mask in integer (0 or 1 value)
        """
        data = np.zeros(width * height, np.uint32)
        im = self.mask_group.get(masktag)
        im2 = snappy.jpy.cast(im, snappy.Mask)  # change from ProductNode to Mask typo
        im2.readPixels(x, y, width, height, data)
        # No need to reshape data as unc values are not!!!
        # data.shape = rectangle.height, rectangle.width
        return data
