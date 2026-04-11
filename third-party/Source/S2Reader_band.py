#!/bin/python3
import os, sys
from glob import glob as glb
import json
import argparse
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate
import rasterio as rio
import xmltodict
from pyproj import Proj
from S2RUT import S2RUT_L1
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ----- S2 reader ------------
S2nodata = 0
S2nodataMax = 10
ndetect = 12

S2bandNames = [
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
S2bandID = {
    "B01": 0,
    "B02": 1,
    "B03": 2,
    "B04": 3,
    "B05": 4,
    "B06": 5,
    "B07": 6,
    "B08": 7,
    "B8A": 8,
    "B09": 9,
    "B10": 10,
    "B11": 11,
    "B12": 12,
}
S2_BAND_SAMPLING = {
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
S2_BAND_PIXELS = {
    "B01": 1830,
    "B02": 10980,
    "B03": 10980,
    "B04": 10980,
    "B05": 5490,
    "B06": 5490,
    "B07": 5490,
    "B08": 10980,
    "B8A": 5490,
    "B09": 1830,
    "B10": 1830,
    "B11": 5490,
    "B12": 5490,
}


class S2Processor:

    def __init__(self):

        # INPUTS -------------------------------------------------------------------------------------------------------
        self.product = None  # path to S2 L1C or S2 L2A product
        self.selected_bands = None  # array with the tags for the bands to be processed
        self.lat_centre = None  # ROI centre. The geographical latitude in decimal degree, valid range is -90 to +90.
        self.lon_centre = None  # ROI centre. The geographical longitude in decimal degree, valid range is -180 to +180.
        self.sensing_time = None  # Timestamp at the S2 L1C acquisition
        self.w = None  # width of the desired data array (meters)
        self.h = None  # height of the desired data array (meters)
        self.noise_model = None  # path to external noise model

        # OUTPUTS ------------------------------------------------------------------------------------------------------
        self.L1C_ref = (
            []
        )  # list of arrays with ROI pixels per band. They provide TOA reflectance values as float.
        self.L1C_rad = []
        self.ds_profile = []
        self.view_azimuth = []  # same for each angle
        self.view_zenith = []
        self.sun_azimuth = []
        self.sun_zenith = []
        self.cloud = None
        self.cirrus = None
        self.metadatadict = {
            "spacecraft": None,
            "quant": None,
            "A": [],
            "offset": [],
            "alpha": {},
            "beta": {},
            "Esun": [],
            "Usun": None,
        }  # this dictionary contains the relevant metadata parameters

    def get_profile(self, product, band):
        try:
            bandfile = glb(f"{product}/GRANULE/*/IMG_DATA/*{band}.jp2")[0]
        except:
            raise ValueError("band ", band, " not found in product", product)

        with rio.open(bandfile) as ds:
            im = ds.read(1).astype(float)
            im[im < S2nodataMax] = np.nan
            return ds.profile

    def read_band(self, product, band):
        try:
            bandfile = glb(f"{product}/GRANULE/*/IMG_DATA/*{band}.jp2")[0]
        except:
            raise ValueError("band ", band, " not found in product", product)

        with rio.open(bandfile) as ds:
            im = ds.read(1).astype(float)
            im[im < S2nodataMax] = np.nan
            return (im - self.metadatadict["offset"][band]) / self.metadatadict[
                "quant"
            ], ds.profile

    def get_win_transform(self, product, band, win):
        try:
            bandfile = glb(f"{product}/GRANULE/*/IMG_DATA/*{band}.jp2")[0]
        except:
            raise ValueError("band ", band, " not found in product", product)

        with rio.open(bandfile) as ds:
            win_transform = ds.window_transform(win)
            # im = ds.read(1, window=win).astype(float)
            # im[im < S2nodataMax] = np.nan
            # im2 = (im - self.metadatadict['offset'][band]) / self.metadatadict['quant']
        return win_transform

    def read_classi(self, product):
        try:
            bandfile = glb(f"{product}/GRANULE/*/QI_DATA/MSK_CLASSI*.jp2")[0]
        except:
            raise ValueError("MSK_CLASSI not found in product", product)

        with rio.open(bandfile) as ds:
            im = ds.read().astype(float)
        return im

    def read_detfoo(self, product, band):
        try:
            bandfile = glb(f"{product}/GRANULE/*/QI_DATA/MSK_DETFOO*{band}.jp2")[0]
        except:
            raise ValueError("MSK_CLASSI not found in product", product)

        with rio.open(bandfile) as ds:
            im = ds.read().astype(float)
        return im

    def read_mtd(self, product, mtd_type):
        try:
            if mtd_type == "MSIL1C":
                mtdfile = glb(product + "/MTD_MSIL1C.xml")[0]
            elif mtd_type == "TL":
                mtdfile = glb(product + "/GRANULE/*/MTD_TL.xml")[0]
            elif mtd_type == "DS":
                mtdfile = glb(product + "/DATASTRIP/*/MTD_DS.xml")[0]
            else:
                raise ValueError(mtd_type, " is not a valid MTD file type")
        except:
            raise ValueError(
                "Could not find MTD file", mtd_type, "for product", product
            )

        with open(mtdfile) as f:
            dic = xmltodict.parse(f.read())
        return dic

    def write_as_band(self, im, product, band, ext):

        # quantize image
        im2 = (im * self.metadatadict["quant"]).astype(int) + self.metadatadict[
            "offset"
        ]
        im2[im == np.nan] = S2nodata

        # get dataset profile
        try:
            bandfile = glb(f"{product}/GRANULE/*/IMG_DATA/*{band}.jp2")[0]
        except:
            raise ValueError("band ", band, " not found in product", product)

        with rio.open(bandfile) as src_ds:
            kwds = src_ds.profile

        kwds.update({"QUALITY": "100", "REVERSIBLE": "YES", "nodata": S2nodata})
        newfile = os.path.basename(bandfile)[:-4] + ext + ".jp2"
        print("Writing ", newfile)
        with rio.open(newfile, "w", **kwds) as dst_ds:
            dst_ds.write(im2, 1)

    def write_unc(self, im, fpath, band, ds_profile):

        # quantize image
        im2 = (
            (im * self.metadatadict["quant"]) + self.metadatadict["offset"][band]
        ).astype("u2")
        im2[im == np.nan] = S2nodata

        # get product path
        # fpath = f'{product}/UNC/'
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        newfile = fpath + band + "_unc.jp2"

        ds_profile["width"] = im2.shape[1]
        ds_profile["height"] = im2.shape[0]

        with rio.open(newfile, "w", **ds_profile) as dst_ds:
            dst_ds.write(im2, 1)

        print("Writing ", newfile)

        # transform = Affine.translation(
        #                         x - res / 2, y - res / 2) * Affine.scale(res, -res)
        #                 with rasterio.open(
        #                 'test.tif',
        #                 'w',
        #                         count=1,
        #                         dtype=np.uint8,
        #                 driver='GTiff',
        #                         height=256,
        #                         width=256,
        #                 crs=rasterio.crs.CRS.from_epsg(3857),
        #                 transform=transform,
        #                 ) as dst:
        #                         dst.write(thearray,1)

    def extrapole_2D_nearest(self, data):
        x = np.arange(0, data.shape[1])
        y = np.arange(0, data.shape[0])
        # mask invalid values
        array = np.ma.masked_invalid(data)
        xx, yy = np.meshgrid(x, y)
        # get only the valid values
        x1 = xx[~array.mask]
        y1 = yy[~array.mask]
        newarr = array[~array.mask]
        GD1 = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method="nearest")
        return GD1

    def get_row_col(self, profile, lat_centre, lon_centre):
        """
        Parameters
        ----------
        profile : band dataset profile, containing the CRS and the transform
        lat_centre : float, decimal latitude
        lon_centre : float, decimal longitude

        Returns
        -------
        ilat : latitude index in band coordinate
        ilon : latitude index in band coordinate
        """
        proj2utm = Proj(profile["crs"], ellps="WGS84")
        utmx, utmy = proj2utm(lon_centre, lat_centre)
        ilat, ilon = rio.transform.rowcol(profile["transform"], utmx, utmy)
        return ilat, ilon

    def get_lat_lon(self, product, band):
        """
        Reconstruct the lat/lon grid of each pixel of a given band from the tile geocoding

        Parameters
        ----------
        product : string, path to S2 L1C product
        band : string, code name of the band.

        Returns
        -------
        lat : 2D array, latitude coordinate in band resolution grid
        lon : 2D array, longitute coordinate in  band resolution grid
        """

        mtd = self.read_mtd(product, "TL")
        epsg = mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
            "HORIZONTAL_CS_CODE"
        ]
        # find the index matching the current band resolution
        for i, l in enumerate(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
                "Geoposition"
            ]
        ):
            if S2_BAND_SAMPLING[band] == int(l["@resolution"]):
                ind = i
        nrows = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"]["Size"][
                ind
            ]["NROWS"]
        )
        ncols = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"]["Size"][
                ind
            ]["NCOLS"]
        )
        ulx = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
                "Geoposition"
            ][ind]["ULX"]
        )
        uly = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
                "Geoposition"
            ][ind]["ULY"]
        )
        xdim = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
                "Geoposition"
            ][ind]["XDIM"]
        )
        ydim = int(
            mtd["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Geocoding"][
                "Geoposition"
            ][ind]["YDIM"]
        )

        x_array = np.linspace(0, nrows, nrows, endpoint=False)
        y_array = np.linspace(0, ncols, ncols, endpoint=False)
        x_grid, y_grid = np.meshgrid(x_array, y_array)
        utm_x = np.add(ulx, np.multiply(x_grid, xdim))
        utm_y = np.add(uly, np.multiply(y_grid, ydim))

        p2utm = Proj(epsg, preserve_units=False)
        lon, lat = p2utm(utm_x, utm_y, inverse=True)

        return lat, lon

    def load_noise_model(self, noise_model, sat):
        """
        Load the alpha and beta coefficient from the specified file instead or
        reading them from the product metadata

        Parameters
        ----------
        noise_model : string, path to noise model file

        """

        # read noise model parameters
        with open(noise_model) as f:
            model = json.load(f)

        for b in self.selected_bands:
            alpha, beta = model[sat]["noise_model"][b]
            # Ak = model['Gains'][b]
            self.metadatadict["alpha"][b] = alpha
            self.metadatadict["beta"][b] = beta
        return

    def get_data(self, theband, get_va=False, doplot=False):
        """
        Parameters
        ----------
        get_va : boolean, optional.
            Activate the computation of viewing angle. The default is False.
        doplot : boolean, optional.
            Activate the plotting of extracted data. The default is False.

        Returns
        -------
        None. (Fill the S2Processor class object with extracted data.)

        """

        ## Get SZA & SAA
        mtd_tl = self.read_mtd(self.product, "TL")
        SZA_LR = np.zeros((23, 23))
        SAA_LR = np.zeros((23, 23))
        SZA_lines = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Angles"][
            "Sun_Angles_Grid"
        ]["Zenith"]["Values_List"]["VALUES"]
        SAA_lines = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"]["Tile_Angles"][
            "Sun_Angles_Grid"
        ]["Azimuth"]["Values_List"]["VALUES"]
        for i in range(len(SZA_lines)):
            SZA_LR[i, :] = np.fromstring(SZA_lines[i], sep=" ")
            SAA_LR[i, :] = np.fromstring(SAA_lines[i], sep=" ")

        # interpolation to HR
        one_axis_LR = np.arange(0, 23 * 5000, 5000)
        SZA_int = interpolate.RectBivariateSpline(one_axis_LR, one_axis_LR, SZA_LR)
        SAA_int = interpolate.RectBivariateSpline(one_axis_LR, one_axis_LR, SAA_LR)

        ## read band profile & detfoo
        rho, profile = self.read_band(self.product, theband)
        # profile = self.get_profile(self.product, theband)
        N = profile["width"]
        res = int(profile["transform"][0])
        detfoo = self.read_detfoo(self.product, theband)[0]

        ## interpolate cos SZA & SAA
        one_axis = np.arange(res / 2, N * res, res)
        SZA_HR = SZA_int(one_axis, one_axis)
        SAA_HR = SAA_int(one_axis, one_axis)

        ## Get VZA and VAA
        if get_va == True:
            VZA_LR = np.zeros((23, 23))
            VAA_LR = np.zeros((ndetect, 23, 23))
            view_angle_list = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"][
                "Tile_Angles"
            ]["Viewing_Incidence_Angles_Grids"]

            # read all data in the list
            for i in range(len(view_angle_list)):
                band_id = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"][
                    "Tile_Angles"
                ]["Viewing_Incidence_Angles_Grids"][i]["@bandId"]
                det_id = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"][
                    "Tile_Angles"
                ]["Viewing_Incidence_Angles_Grids"][i]["@detectorId"]

                # if band_id matchs the current band: read angles
                if int(band_id) == S2bandID[theband]:
                    VZA = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"][
                        "Tile_Angles"
                    ]["Viewing_Incidence_Angles_Grids"][i]["Zenith"]["Values_List"][
                        "VALUES"
                    ]
                    VAA = mtd_tl["n1:Level-1C_Tile_ID"]["n1:Geometric_Info"][
                        "Tile_Angles"
                    ]["Viewing_Incidence_Angles_Grids"][i]["Azimuth"]["Values_List"][
                        "VALUES"
                    ]
                    idet = int(det_id) - 1
                    for j in range(len(VZA)):
                        # VZA: keep finite value from each detector (0 otherwise)
                        vza = np.fromstring(VZA[j], sep=" ")
                        VZA_LR[j, :] = np.where(np.isfinite(vza), vza, VZA_LR[j, :])
                        # VAA: store data according to the detector index
                        VAA_LR[idet, j, :] = np.fromstring(VAA[j], sep=" ")

                    # VAA: replace nan with basic extrapolation for viewing azimuth angle.
                    VAA_LR[idet, :, :] = self.extrapole_2D_nearest(VAA_LR[idet, :, :])
                    # VAA_LR[idet,:,:] = np.where(np.isfinite(VAA_LR[idet,:,:]), VAA_LR[idet,:,:], np.nanmedian(VAA_LR[idet,:,:]))

            # VZA: Interpolate to HR directly
            VZA_HR = np.zeros((N, N))
            VZA_int = interpolate.RectBivariateSpline(one_axis_LR, one_axis_LR, VZA_LR)
            VZA_HR = VZA_int(one_axis, one_axis)
            # VAA: loop on detector to properly interpolate them at HR by taking into account det footprint
            VAA_HR = np.zeros((N, N))
            for i in range(ndetect):
                if i + 1 in detfoo:
                    # VAA: interpolate single detector data to HR and keep only the actual footprint
                    VAA_int = interpolate.RectBivariateSpline(
                        one_axis_LR, one_axis_LR, VAA_LR[i, :, :]
                    )
                    VAA_HR_single_det = VAA_int(
                        one_axis, one_axis
                    )  # single detector data fully interpolated
                    # return data from single detector if condition is true, otherwise return the original data
                    VAA_HR = np.where(detfoo == i + 1, VAA_HR_single_det, VAA_HR)

            if doplot:
                plt.figure()
                plt.imshow(VZA_HR)
                plt.title("VZA " + theband)
                plt.colorbar(label="Angle (°)")
                plt.figure()
                plt.imshow(VAA_HR)
                plt.title("VAA " + theband)
                plt.colorbar(label="Angle (°)")

        ## ROI centre: If no lat/lon specified set ROI in tile centre, else get row/col index of ROI centre
        if (self.lat_centre is None) or (self.lon_centre) is None:
            ilat, ilon = N / 2, N / 2
        else:
            ilat, ilon = self.get_row_col(profile, self.lat_centre, self.lon_centre)

        ## ROI size: if none return whole image.
        if (self.w is None) or (self.h is None):
            xstart, xstop = 0, N - 1
            ystart, ystop = 0, N - 1
            npix_x, npix_y = N / 2, N / 2
        else:
            npix_x = self.h // res
            xstart, xstop = ilat - npix_x / 2, ilat + npix_x / 2
            npix_y = self.w // res
            ystart, ystop = ilon - npix_y / 2, ilon + npix_y / 2

            # if roi size exceed product boudaries return whole image.
            if xstart < 0 or xstop >= N:
                print("WARNING: ROI x-axis size exceed product boundaries !!")
                xstart, xstop = 0, N - 1
            if ystart < 0 or ystop >= N:
                print("WARNING: ROI y-axis size exceed product boundaries !!")
                ystart, ystop = 0, N - 1

        ## Read 'transform' function of selected windows for geoloc
        xstart, xstop = int(xstart), int(xstop)
        ystart, ystop = int(ystart), int(ystop)
        win = rio.windows.Window(ystart, xstart, npix_y, npix_x)
        win_trans = self.get_win_transform(self.product, theband, win)

        ## update profile with ROI 'transform'
        profile["transform"] = win_trans

        ## compute radiance L
        L = (
            rho[xstart:xstop, ystart:ystop]
            * self.metadatadict["Esun"][theband]
            * self.metadatadict["Usun"]
            * np.cos(SZA_HR[xstart:xstop, ystart:ystop] * np.pi / 180.0)
            / np.pi
        )

        if doplot:
            plt.figure()
            plt.imshow(rho)
            plt.plot(ilon, ilat, "*r")
            plt.title(theband + " reflectance")
            plt.colorbar(label="Reflectance")

            plt.figure()
            plt.imshow(L)
            plt.title(theband + " ROI radiance")
            plt.colorbar(label="Radiance")

        # # output value
        # self.L1C_ref.append(rho[xstart:xstop,ystart:ystop])
        # self.L1C_rad.append(L)
        self.sun_azimuth.append(SAA_HR[xstart:xstop, ystart:ystop])
        self.sun_zenith.append(SZA_HR[xstart:xstop, ystart:ystop])
        self.ds_profile.append(profile)
        # if get_va:
        #     self.view_azimuth.append(VAA_HR[xstart:xstop,ystart:ystop])
        #     self.view_zenith.append(VZA_HR[xstart:xstop,ystart:ystop])

        # clear variable for memory usage
        del rho
        del SAA_HR, SAA_int
        del SZA_HR, SZA_int
        if get_va:
            del VAA_HR, VAA_int, VAA_LR
            del VZA_HR, VZA_int, VZA_LR

        return L, profile

    def get_metadata(self):

        sat = os.path.basename(self.product)[0:3]

        # get E0, U, quant & offset
        mtd_eup = self.read_mtd(self.product, "MSIL1C")
        img_data = mtd_eup["n1:Level-1C_User_Product"]["n1:General_Info"][
            "Product_Image_Characteristics"
        ]
        quant = float(img_data["QUANTIFICATION_VALUE"]["#text"])

        offset = {}
        for x in img_data["Radiometric_Offset_List"]["RADIO_ADD_OFFSET"]:
            band = S2bandNames[int(x["@band_id"])]
            if band in self.selected_bands:
                offset[band] = float(x["#text"])

        refl_data = img_data["Reflectance_Conversion"]
        U = float(refl_data["U"])
        E0 = {}
        for x in refl_data["Solar_Irradiance_List"]["SOLAR_IRRADIANCE"]:
            band = S2bandNames[int(x["@bandId"])]
            if band in self.selected_bands:
                E0[band] = float(x["#text"])

        # get physical gain
        mtd_ds = self.read_mtd(self.product, "DS")
        spect_data = mtd_ds["n1:Level-1C_DataStrip_ID"]["n1:Image_Data_Info"][
            "Sensor_Configuration"
        ]["Acquisition_Configuration"]["Spectral_Band_Info"][
            "Spectral_Band_Information"
        ]
        A_gain = {}
        for x in spect_data:
            band = S2bandNames[int(x["@bandId"])]
            if band in self.selected_bands:
                A_gain[band] = float(x["PHYSICAL_GAINS"]["#text"])

        # Load noise model
        if self.noise_model is not None:
            self.load_noise_model(self.noise_model, sat)
        else:
            # if not defined, read default noise model but warn users.
            print("WARNING: Noise model is outdated and not representative !")
            print("         You should use an extrernal noise model instead.")
            rad_data = mtd_ds["n1:Level-1C_DataStrip_ID"]["n1:Quality_Indicators_Info"][
                "Radiometric_Info"
            ]["Radiometric_Quality_List"]["Radiometric_Quality"]
            alpha_dict = {}
            beta_dict = {}
            for x in rad_data:
                band = S2bandNames[int(x["@bandId"])]
                if band in self.selected_bands:
                    alpha_dict[band] = float(x["Noise_Model"]["ALPHA"])
                    beta_dict[band] = float(x["Noise_Model"]["BETA"])

            self.metadatadict["alpha"] = alpha_dict
            self.metadatadict["beta"] = beta_dict

        # store metadata
        if os.path.basename(self.product)[0:3] == "S2A":
            self.metadatadict["spacecraft"] = "Sentinel-2A"
        elif os.path.basename(self.product)[0:3] == "S2B":
            self.metadatadict["spacecraft"] = "Sentinel-2B"
        elif os.path.basename(self.product)[0:3] == "S2C":
            self.metadatadict["spacecraft"] = "Sentinel-2C"

        self.metadatadict["quant"] = quant
        self.metadatadict["Usun"] = U
        self.metadatadict["Esun"] = E0
        self.metadatadict["A"] = A_gain
        self.metadatadict["offset"] = offset

        # get mask classification
        msk_classi = self.read_classi(self.product)
        self.cloud = msk_classi[0, :, :]
        self.cirrus = msk_classi[1, :, :]
        # msk_snow = msk_classi[2,:,:]
        del msk_classi

    def get_unc(self, get_va=False, doplot=False):

        RUTl1 = S2RUT_L1()

        ## loop over bands
        for ib, theband in enumerate(self.selected_bands):
            print("Extracting:", theband)
            L, prof = self.get_data(theband, get_va, doplot)
            u_ref, u_sig, u_sys, u_cont = RUTl1.unc_calculation(
                L, theband, ib, self.metadatadict, self.sun_zenith
            )

            # write unc product
            self.write_unc(u_ref, self.output, theband, self.ds_profile[ib])

            # Plot
            if doplot:
                fig1, ax1 = plt.subplots(1, 2, figsize=(8, 6))
                im0 = ax1[0].imshow(L, interpolation="None")
                ax1[0].set_title(theband + " Radiance")
                divider = make_axes_locatable(ax1[0])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig1.colorbar(im0, cax=cax, orientation="vertical", label="radiance")

                im1 = ax1[1].imshow(u_ref, interpolation="None")
                ax1[1].set_title(theband + " Uncertainties")
                divider = make_axes_locatable(ax1[1])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig1.colorbar(im1, cax=cax, orientation="vertical", label="unc [%]")
