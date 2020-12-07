# -*- coding: utf-8 -*-
"""
Almost working to put into Java but not quite - fails on the setPixels function

EOpy adaptation of the Sentinel-2 Radiometric Uncertainty tool

Based on the python code found in the git repository:
https://github.com/senbox-org/snap-rut/tree/master/src/main/python

- I've not changed any of the calculations
- I've only made it eopy compatible

Created 04/03/2020 by N. Origo (original code by Javier Gorrono)
"""
import snappy
from snappy import HashMap as hash
import s2_rut_python.s2_rut_algo
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

        # make params visible to the other functions
        self.params = params
        self.product = product

        self.source_product = product.product[0]["product"]
        """
        # compute for all resolutions
        for p in range(3):
            self.source_product   = product.product[p]["product"]
            self.initialise()

            for tb in self.toa_band_names:
                #band        = self.source_product.getBand(tb)
                bandu       = self.rut_product.getBand(tb+"_rut")
                #srci        = band.getSourceImage()
                #tiles       = [srci.getTiles()]
                uband,bw,bh = self.computeBand(bandu)
                ubandrs     = uband.reshape(bh,bw)
                # add the band to the in memory rut product
                print p,tb
                bandu.ensureRasterData()
                for y in range(bh):
                    bandu.writePixels(0,y,bw,1,uband[y*bw:y*bw+bw])
                #bandu.setPixels(0,0,bw,bh,uband)

                # maybe we don't even need to do this
                #for ti,tile in enumerate(tiles):
                #    # here we need to put the tile into the java product
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
        granules_meta = metadata_root.getElement("Granules")

        self.spacecraft = (
            self.datastrip_meta.getElement("General_Info")
            .getElement("Datatake_Info")
            .getAttributeString("SPACECRAFT_NAME")
        )

        # todo - check if there is a granule

        if "band_names" in self.params.keys():
            self.toa_band_names = self.params["band_names"]
        else:
            self.toa_band_names = S2_BAND_NAMES

        self.rut_algo.u_sun = self.get_u_sun(self.product_meta)
        self.rut_algo.quant = self.get_quant(self.product_meta)

        self.source_sza = self.source_product.getBand("sun_zenith")

        self.rut_algo.k = self.get_k()
        self.rut_algo.unc_select = self.get_unc_select()

        self.sourceBandMap = {}
        for name in self.toa_band_names:
            # TODO - Change the interface so that undesired bands (e.g azimuth) are not shown.
            if (
                not name in S2_BAND_NAMES
            ):  # The band name is checked to confirm it is valid band.
                if ("view_" in name) or ("sun_" in name):
                    continue  # the angular bands are shown in the GUI and we simply jump to the next band if selected
                else:
                    raise RuntimeError(
                        'Source band "'
                        + name
                        + '" is not valid and has not been processed'
                    )

            source_band = self.source_product.getBand(name)
            unc_toa_band = snappy.Band(
                name + "_rut",
                snappy.ProductData.TYPE_UINT8,
                source_band.getRasterWidth(),
                source_band.getRasterHeight(),
            )
            unc_toa_band.setDescription(
                "Uncertainty of "
                + name
                + " (coverage factor k="
                + str(self.rut_algo.k)
                + ")"
            )
            unc_toa_band.setNoDataValue(250)
            unc_toa_band.setNoDataValueUsed(True)
            self.targetBandList.append(unc_toa_band)
            self.sourceBandMap[unc_toa_band] = source_band
            snappy.ProductUtils.copyGeoCoding(source_band, unc_toa_band)

        masterband = self.get_masterband(self.targetBandList)
        rut_product = snappy.Product(
            self.source_product.getName() + "_rut",
            "S2_RUT",
            masterband.getRasterWidth(),
            masterband.getRasterHeight(),
        )  # in-memory product
        snappy.ProductUtils.copyGeoCoding(masterband, rut_product)
        for band in self.targetBandList:
            rut_product.addBand(band)

        # The metadata from the RUT product is defined
        self.rut_product_meta = (
            rut_product.getMetadataRoot()
        )  # Here we define the product metadata
        # SOURCE_PRODUCT
        sourceelem = MetadataElement("Source_product")
        data = snappy.ProductData.createInstance(self.source_product.getDisplayName())
        sourceattr = MetadataAttribute(
            "SOURCE_PRODUCT", snappy.ProductData.TYPE_ASCII, data.getNumElems()
        )
        sourceattr.setData(data)
        sourceelem.addAttribute(sourceattr)
        self.rut_product_meta.addElement(sourceelem)
        # COVERAGE_FACTOR
        sourceelem = MetadataElement("Coverage_factor")
        data = snappy.ProductData.createInstance(str(self.rut_algo.k))
        sourceattr = MetadataAttribute(
            "COVERAGE_FACTOR", snappy.ProductData.TYPE_ASCII, data.getNumElems()
        )
        sourceattr.setData(data)
        sourceelem.addAttribute(sourceattr)
        self.rut_product_meta.addElement(sourceelem)
        # RUT_VERSION
        version = (
            snappy.GPF.getDefaultInstance()
            .getOperatorSpiRegistry()
            .getOperatorSpi("S2RutOp")
            .getOperatorDescriptor()
            .getVersion()
        )
        sourceelem = MetadataElement("Version")
        data = snappy.ProductData.createInstance(version)
        sourceattr = MetadataAttribute(
            "VERSION", snappy.ProductData.TYPE_ASCII, data.getNumElems()
        )
        sourceattr.setData(data)
        sourceelem.addAttribute(sourceattr)
        self.rut_product_meta.addElement(sourceelem)
        # CONTRIBUTOR LIST: List of selected ones
        sourceelem = MetadataElement("List_Contributors")
        contributors = [
            "INSTRUMENT_NOISE",
            "OOF_STRAYLIGHT-SYSTEMATIC",
            "OOF_STRAYLIGHT-RANDOM",
            "CROSSTALK",
            "ADC_QUANTISATION",
            "DS_STABILITY",
            "GAMMA_KNOWLEDGE",
            "DIFFUSER-ABSOLUTE_KNOWLEDGE",
            "DIFFUSER-TEMPORAL_KNOWLEDGE",
            "DIFFUSER-COSINE_EFFECT",
            "DIFFUSER-STRAYLIGHT_RESIDUAL",
            "L1C_IMAGE_QUANTISATION",
        ]
        for i in range(0, len(contributors)):
            data = snappy.ProductData.createInstance(str(self.rut_algo.unc_select[i]))
            sourceattr = MetadataAttribute(
                contributors[i], snappy.ProductData.TYPE_ASCII, data.getNumElems()
            )
            sourceattr.setData(data)
            sourceelem.addAttribute(sourceattr)
        self.rut_product_meta.addElement(sourceelem)
        # DATE OF PROCESSING
        sourceelem = MetadataElement("Processing_datetime")
        data = snappy.ProductData.createInstance(str(datetime.datetime.now()))
        sourceattr = MetadataAttribute(
            "PROCESSING_DATETIME", snappy.ProductData.TYPE_ASCII, data.getNumElems()
        )
        sourceattr.setData(data)
        sourceelem.addAttribute(sourceattr)
        self.rut_product_meta.addElement(sourceelem)

        # context.setTargetProduct(rut_product)
        self.rut_product = rut_product

    def computeBand(self, bandu):
        # Logging template
        # SystemUtils.LOG.info('target band name: ' + band.getName())
        # SystemUtils.LOG.info('tile rect: ' + tile.getRectangle().toString())
        source_band = self.sourceBandMap[bandu]
        gname = source_band.getName()
        toa_band_id = np.int(S2_BAND_NAMES.index(gname))

        source_sza = self.source_sza

        bw = source_band.getRasterWidth()
        bh = source_band.getRasterHeight()

        for resolution in np.arange(3):
            pprod = self.product.product[resolution]
            if pprod["product"] == source_band.getProduct():
                cloudmask = self.mask_roi(
                    "opaque_clouds_" + pprod["product_name"].split("_")[0], 0, 0, bw, bh
                )
                cirrusmask = self.mask_roi(
                    "cirrus_clouds_" + pprod["product_name"].split("_")[0], 0, 0, bw, bh
                )
                break

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

        degrademask = self.mask_roi("msi_degraded_" + gname, 0, 0, bw, bh)
        lostmask = self.mask_roi("msi_lost_" + gname, 0, 0, bw, bh)
        defectmask = self.mask_roi("defective_" + gname, 0, 0, bw, bh)
        invalidmask = np.maximum(np.maximum(degrademask, lostmask), defectmask)
        satl1amask = self.mask_roi("saturated_l1a_" + gname, 0, 0, bw, bh)
        satl1bmask = self.mask_roi("saturated_l1b_" + gname, 0, 0, bw, bh)
        satl1mask = np.maximum(satl1amask, satl1bmask)
        nodatamask = self.mask_roi("nodata_" + gname, 0, 0, bw, bh)
        # selects the maximum element-wise. Mask true value is 255. This is higher than 250 (max uncertainty permitted)
        # 251 is for degraded,lost or defective data. 252 is for saturated (L1a or L1b). 253 is for pixel with no data,
        # 254 is for cirrus cloud and 255 is for opaque clouds.
        val = np.maximum(unc, np.uint8(251 * invalidmask / 255))
        val = np.maximum(val, np.uint8(252 * satl1mask / 255))
        val = np.maximum(val, np.uint8(253 * nodatamask / 255))
        val = np.maximum(val, np.uint8(254 * cirrusmask / 255))
        val = np.maximum(val, np.uint8(cloudmask))
        # tile.setSamples(val)
        return val, bw, bh

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
