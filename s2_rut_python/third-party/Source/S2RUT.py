# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:48:33 2016
@author: jg9
"""

import math
import numpy as np
from scipy.stats import multivariate_normal
import json
from matplotlib import pyplot as plt


# S2_BAND_NAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
S2_BAND_NAMES = ["B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B10","B11","B12"]
S2_BAND_SAMPLING = {'B01':60,'B02':10,'B03':10,'B04':10,'B05':20,'B06':20,'B07':20,'B08':10,'B8A':20,'B09':60,'B10':60,'B11':20,'B12':20}

class S2RUT_L1:
    
    def __init__(self, contributors_file):
        """
        Contains Sentinel-2 Level-1 radiometric configuration and uncertainty contributions placeholders.
        """
        with open(contributors_file) as f:
            contributors = json.load(f)

        # read contributors value
        self.Lref = contributors["Lref"]
        self.u_stray_rand_all = contributors["u_stray_rand_all"]
        self.u_xtalk_all = contributors["u_xtalk_all"] # units in W.m-2.sr-1.μm-1
        self.u_DS_all = contributors["u_DS_all"]
        self.u_diff_absarray = contributors["u_diff_absarray"] # % values from ICCDB (S2A at S2_OPER_MSI_DIFF___20150519T000000_0009.xml and S2B at S2_OPER_MSI_DIFF___20160415T000000_0001.xml)
        self.u_diff_cos = contributors["u_diff_cos"] # [%] from 0.13° diffuser planarity/micro as in (AIRBUS 2015). Assumed same for S2A/S2B.
        self.u_diff_sl = contributors["u_diff_sl"] # [%] as a conservative residual (AIRBUS 2015). Assumed same for S2A/S2B.
        self.u_diff_temp = contributors["u_diff_temp"]  # This value is correctly redefined for specific satellite at the S2RutOp.
        self.u_ADC = contributors["u_ADC"] # [DN](rectangular distribution, see combination)
        self.u_gamma = contributors["u_gamma"]
        self.k = contributors["k"]  # This value is correctly redefined for specific satellite at the S2RutOp.
        self.u_resamp = contributors["u_resamp"]  # noise resampling factor from L1b to L1c (https://doi.org/10.1117/12.2603730)
        # self.u_geoloc = contributors["u_geoloc"]  # geoloc err automatically defined. 1.5m or 3m depending on refinement
        
        # list of booleans with user selected uncertainty sources
        self.unc_select_noise = contributors["unc_select_noise"]
        self.unc_select_stray_sys = contributors["unc_select_stray_sys"]
        self.unc_select_stray_rand = contributors["unc_select_stray_rand"]
        self.unc_select_xtalk = contributors["unc_select_xtalk"]
        self.unc_select_adc = contributors["unc_select_adc"]
        self.unc_select_ds = contributors["unc_select_ds"]
        self.unc_select_gamma = contributors["unc_select_gamma"]
        self.unc_select_diff_abs = contributors["unc_select_diff_abs"]
        self.unc_select_diff_temp = contributors["unc_select_diff_temp"]
        self.unc_select_diff_cos = contributors["unc_select_diff_cos"]
        self.unc_select_diff_sl = contributors["unc_select_diff_sl"]
        self.unc_select_ref_quant = contributors["unc_select_ref_quant"]
        self.unc_select_geoloc= contributors["unc_select_geoloc"]
        


    def unc_calculation_abs(self, band_data, bandname, bandind, metadatadict, sun_zenith, do_contributor=False):
        """
        This function represents the core of the RUTv1. It takes as an input the pixel data of a specific band and tile
        in a S2-L1C product and produces an image with the same dimensions that contains the radiometric uncertainty of
        each pixel reflectance factor. The steps and its numbering is equivalent to the RUT-DPM. This document can be
        found in the tool github. Also there a more detailed explanation of the theoretical background can be found.
        [aderu 2025-02-24] The function now work in absolute unit
        :param band_data: list with the quantized L1C reflectance pixels of a band
        :param bandname: tag name of the S2 L1C band
        :param bandind: zero-based index of the band in the current processed list
        :param metadatadict: dictionary contains the relevant metadata parameters
        :param sun_zenith: sunzenith angle per pixel level for each band
        :param do_contributor: activate the per contributor output 
        
        :return: u_ref: uncertainty associated to each pixel in reflectance unit.
        :return: u_contributor: struct with uncertainty for each contributor in reflectance unit.
        """
        
        band_id = np.where(np.array(S2_BAND_NAMES) == bandname)[0][0]  # check for index in full L1C band names!
       
        #######################################################################
        # 0.	set Negative radiance to 0
        #######################################################################
        band_data[band_data<0] = 0
        
        #######################################################################
        # 1.	Undo reflectance conversion
        #######################################################################
        # a.	No action required
        # b.	[product metadata] #issue: missing one band
        #    General_Info/Product_Image_Characteristics/PHYSICAL_GAINS [bandId]
        #    [datastrip metadata]
        #    Image_Data_Info/Sensor_Configuration/Acquisition_Configuration/
        #    Spectral_Band_Info/Spectral_Band_Information [bandId]/ PHYSICAL_GAINS

        # Replace the reflectance factors by CN values
        cn = (metadatadict['A'][bandname] * metadatadict['Esun'][bandname] * metadatadict['Usun'] * np.cos(np.radians(sun_zenith[bandind])) / math.pi) * band_data

        #######################################################################
        # 2.	Orthorectification process
        #######################################################################

        # geoloc error (computed from LSB gradient)
        if self.unc_select_geoloc:
            # geo = self.u_geoloc/S2_BAND_SAMPLING[bandname]
            if metadatadict['refined'] == True:
                # 1.5 meters for refined tiles
                geo = 1.5/S2_BAND_SAMPLING[bandname] 
            else:
                # 3 meters for non-refined tiles
                geo = 3/S2_BAND_SAMPLING[bandname]
        else:
            geo = 0 
            
        u_geo = geo*np.sqrt(np.gradient(cn)[0]**2 + np.gradient(cn)[1]**2)
           
        #######################################################################
        # 3.	L1B uncertainty contributors: raw and dark signal
        #######################################################################

        # Noise model and Resampling from L1b to L1c (noise model applied to LSB value)
        if self.unc_select_noise:
            u_noise = self.u_resamp * np.sqrt(metadatadict['alpha'][bandname]**2 + metadatadict['beta'][bandname]*cn)
        else:
            u_noise = 0

        # Systematic out-of-field straylight 0.3%*Lmean (instead of Lref) (AIRBUS 2015 & OPT-MPC 2025)
        if self.unc_select_stray_sys:
            # stray_sys = 0.3 * self.Lref[band_id] / 100
            # u_stray_sys = np.ones(np.shape(cn))*metadatadict['A'][bandname] * stray_sys
            u_stray_sys = 0.3 * np.nanmean(cn) * np.ones(np.shape(cn)) / 100
            u_stray_sys[~np.isfinite(cn)] = np.nan # reproduce nan data from input
        else:
            u_stray_sys = 0

        # random straylight [%] from (AIRBUS 2015) and (AIRBUS 2012), compute absolute LSB
        if self.unc_select_stray_rand:
            u_stray_rand = self.u_stray_rand_all[metadatadict['spacecraft']][band_id]*cn/100  
        else:
            u_stray_rand = 0

        # Cross talk [W.m-2.sr-1.μm-1] (AIRBUS 2015), converted to LSB with physical gain.
        if self.unc_select_xtalk:
            u_xtalk = metadatadict['A'][bandname] * self.u_xtalk_all[metadatadict['spacecraft']][band_id]
        else:
            u_xtalk = 0

        # ADC quantisation
        if self.unc_select_adc:
            u_adc = np.ones(np.shape(cn))*self.u_ADC/math.sqrt(3) 
            u_adc[~np.isfinite(cn)] = np.nan # reproduce nan data from input

        else:
            u_adc = 0  

        # Dark signal [LSB]
        if self.unc_select_ds:
            u_DS = self.u_DS_all[metadatadict['spacecraft']][band_id]
            u_ds = np.ones(np.shape(cn))*u_DS
            u_ds[~np.isfinite(cn)] = np.nan # reproduce nan data from input
        else:
            u_ds = 0

        #######################################################################
        # 4.	L1B uncertainty contributors: gamma correction
        #######################################################################
        if self.unc_select_gamma:
            u_gamma = self.u_gamma[band_id]*cn/100
        else:
            u_gamma = 0  # predefined but updated to 0 if deselected by user

        #######################################################################
        # 5.	L1C uncertainty contributors: absolute calibration coefficient
        #######################################################################
        
        # absolute BRDF accuracy, [%] from ICCDB
        if self.unc_select_diff_abs:
            u_diff_abs = self.u_diff_absarray[metadatadict['spacecraft']][band_id]*cn/100
        else:
            u_diff_abs = 0

        if self.unc_select_diff_temp:
            u_diff_temp = self.u_diff_temp*cn/100
        else:
            u_diff_temp = 0

        # diffuser cosine effect [%] 
        if self.unc_select_diff_cos:
            u_diff_cos = self.u_diff_cos*cn/100
        else:                
            u_diff_cos = 0 
            

        if self.unc_select_diff_sl:
            u_diff_sl = self.u_diff_sl*cn/100
        else:
            u_diff_sl = 0  

        #######################################################################
        # 6.	L1C uncertainty contributors: reflectance conversion
        #######################################################################
        if self.unc_select_ref_quant:
            u_ref_quant = np.ones(np.shape(cn)) * (0.5/math.sqrt(3)) / (metadatadict['quant'])  # [%] scaling 0-1 in steps number=quant
            u_ref_quant[~np.isfinite(cn)] = np.nan # reproduce nan data from input
        else:
            u_ref_quant = 0

        #######################################################################
        # 7.	Combine uncertainty contributors
        #######################################################################
        # u_1sigma & u_sys are intermediate steps not computed to save RAM.

        u_stray = np.sqrt(u_stray_rand**2 + u_xtalk**2)
        u_diff = np.sqrt(u_diff_abs**2 + u_diff_cos**2 + u_diff_sl**2)
        # u_1sigma = np.sqrt(u_ref_quant**2 + u_gamma**2 + u_stray**2 + u_diff**2 + u_noise**2 + u_adc**2 + u_ds**2 + u_geo**2)
        # u_sys = u_diff_temp + u_stray_sys
        
        # absolute computation in CN, factor needed to convert into reflectance dimension
        ref_conv = np.pi/(metadatadict['A'][bandname]*metadatadict['Esun'][bandname]*metadatadict['Usun']*np.cos(np.radians(sun_zenith[bandind])))
        
        # total unc: u = u_sys + k * u_1sigma
        u_ref = ( (u_diff_temp+u_stray_sys) + self.k *
                 np.sqrt(u_ref_quant**2+u_gamma**2+u_stray**2+u_diff**2+u_noise**2+u_adc**2+u_ds**2+u_geo**2)
                 ) * ref_conv
        
       
        # export all contributors (if selected by user)
        u_contributor = {}
        if do_contributor:
            if self.unc_select_noise: u_contributor["u_noise"] = u_noise * ref_conv
            if self.unc_select_stray_sys or self.unc_select_diff_temp: 
                u_contributor["u_sys"] = (u_diff_temp+u_stray_sys) * ref_conv
            if self.unc_select_stray_rand or self.unc_select_xtalk: 
                u_contributor["u_stray_xtalk"] = u_stray * ref_conv
            if self.unc_select_adc: u_contributor["u_adc"] = u_adc * ref_conv
            if self.unc_select_ds: u_contributor["u_ds"] = u_ds * ref_conv
            if self.unc_select_gamma: u_contributor["u_gamma"] = u_gamma * ref_conv
            if self.unc_select_diff_abs or self.unc_select_diff_cos or self.unc_select_diff_sl:
                u_contributor["u_diff"] = u_diff * ref_conv
            if self.unc_select_ref_quant: u_contributor["u_ref_quant"] = u_ref_quant * ref_conv
            if self.unc_select_geoloc: u_contributor["u_geoloc"] = u_geo * ref_conv
                                    
        return u_ref, u_contributor

    

    def unc_spectralcorrelation(self, band_rad, metadatadict, rep):
        """
        :param band_data: list with the quantized L1C reflectance values for each band (e.g. at Lref)
        :param band_rad: list with the quantized L1C radiance values for each band (e.g. Lref)
        :param metadatadict: dictionary contains the relevant metadata parameters for typically S2 L2A bands
        :return: samples out of the potential uncertainty L1C for each band. list of arrays [band, #samples]
        """


        def correlated_samples(means, cov_matrix, numbersamples):
            '''
            E.g. for the case of two variables it would be
            cov_matrix = [[(SIGMA1 ** 2), CORR * SIGMA1 * SIGMA2], [CORR * SIGMA1 * SIGMA2, (SIGMA2 ** 2)]]
            distribution = multivariate_normal(mean=[MU1, MU2], cov=cov_matrix)
            samples = distribution.rvs(1000000)
            samples1 = samples[:, 0]
            sample2 = samples[:, 1]
            # A nice example can be found here.
            https://towardsdatascience.com/correlated-variables-in-monte-carlo-simulations-19266fb1cf29
            :param means: list of means e.g. [MU1, MU2].
            :param cov_matrix: covariance matrix (covariance is equal to CORR * SIGMA1 * SIGMA2)
            :param numbersamples: integer with the number of samples that are required.
            :return: array of samples with shape [numbersamples, number variables]
            '''
            # TODO - a multivariate normal is a first approach but there are many other options that will need to consider.
            # see more here https://docs.scipy.org/doc/scipy/reference/stats.html#multivariate-distributions
            distribution = multivariate_normal(mean=means, cov=cov_matrix, allow_singular=True)
            return distribution.rvs(numbersamples)

        # ADERU: had to change to extract the Acoeff vector for all band dictionnary
        Acoeff = []
        alpha = []
        beta = []
        for k in metadatadict['A'].keys():
            Acoeff.append(metadatadict['A'][k])
            alpha.append(metadatadict['alpha'][k])
            beta.append(metadatadict['beta'][k])

        #######################################################################
        # 1.	Undo reflectance conversion
        #######################################################################
        cn = Acoeff * np.array(band_rad)
        # cn = np.array(metadatadict['A']) * band_rad
        #######################################################################
        # 2.	Orthorectification process
        #######################################################################
        # TBD. Here both terms will be used with no distinction.
        #######################################################################
        # 3.	L1B uncertainty contributors: raw and dark signal
        #######################################################################
        u_expand = []
        # noise = 100 * np.sqrt(np.array(metadatadict['alpha']) ** 2 + np.array(metadatadict['beta']) * cn) / cn
        noise = 100 * np.sqrt(np.array(alpha) ** 2 + np.array(beta) * cn) / cn
        Lref = np.hstack((self.Lref[0:10], self.Lref[11:13]))  # We remove B10
        u_stray_sys = 0.3 * np.array(Lref) / band_rad  # adapted to L2A bands 0.3% of Lref.
        # [W.m-2.sr-1.μm-1] 0.3%*Lref all bands (AIRBUS 2015) and (AIRBUS 2014)
        # This a known bias added linearly. when a e.g. band ratio is set, the term Lref_i/Lref_j should be kept.
        stray_rand = np.hstack((self.u_stray_rand_all[metadatadict['spacecraft']][0:10],
                                self.u_stray_rand_all[metadatadict['spacecraft']][
                                11:13]))  # [%](AIRBUS 2015) and (AIRBUS 2012)
        xtalk = np.hstack((self.u_xtalk_all[metadatadict['spacecraft']][0:10],
                           self.u_xtalk_all[metadatadict['spacecraft']][11:13]))  # [W.m-2.sr-1.μm-1](AIRBUS 2015)
        u_DSunc = np.hstack((self.u_DS_all[metadatadict['spacecraft']][0:10],
                             self.u_DS_all[metadatadict['spacecraft']][11:13]))

        corr_matrix = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B1 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B2 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B3 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B4 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B5 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B6 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B7 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B8 vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B8A vs B1-12
                               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],  # B9 vs B1-12
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],  # B11 vs B1-12
                               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])  # B12 vs B1-12
        U = np.array(
            self.u_diff_absarray[metadatadict['spacecraft']][0:10] + self.u_diff_absarray[metadatadict['spacecraft']][
                                                                     11:13])
        cov_diff_abs = np.dot(U.reshape(-1,1), U.reshape(-1,1).T) * corr_matrix

        for i in range(rep):
            u_noise = np.random.normal(0, noise)  # uncorrelated samples for each band
            u_stray_rand = np.random.normal(0, stray_rand)  # uncorrelated samples for each band
            u_xtalk = np.random.normal(0, xtalk)  # uncorrelated samples for each band
            u_DS = np.random.normal(0, u_DSunc)  # correlated samples for each band
            #######################################################################
            # 4.	L1B uncertainty contributors: gamma correction
            #######################################################################
            u_gamma = correlated_samples(np.zeros(12),  (self.u_gamma**2 * corr_matrix), 1)  # [%] (AIRBUS 2015)
            #######################################################################
            # 5.	L1C uncertainty contributors: absolute calibration coefficient
            #######################################################################
            u_diff_abs = correlated_samples(np.zeros(12), cov_diff_abs, 1)
            #######################################################################
            # 6.	L1C uncertainty contributors: reflectance conversion
            #######################################################################
            # not included as a simplification. Small impact under really low reflectance e.g. 0.01 reflectance, 0.3% unc.
            # u_ref_quant = 100 * (0.5 / math.sqrt(3)) / (metadatadict['quant'] * band_data)  # [%]scaling 0-1 in steps number=quant
            #######################################################################
            # 7.	Combine uncertainty contributors
            #######################################################################
            # NOTE: no gamma propagation for RUTv1!!!
            # values given as percentages. Multiplied by 10 and saved to 1 byte(uint8)
            # Clips values to 0-250 --> uncertainty >=25%  assigns a value 250.
            # Uncertainty <=0 represents a processing error (uncertainty is positive)
            u_adc = 100 * np.random.uniform(-self.u_ADC, self.u_ADC, 12) / cn
            u_ds = (100 * u_DS) / cn

            u_stray = u_stray_rand * ((100 * np.array(Acoeff) * u_xtalk) / cn)
            u_diff = u_diff_abs * np.random.normal(0, self.u_diff_cos, 12) + np.random.normal(0, self.u_diff_sl, 12)
            u_1sigma = u_gamma + u_stray + u_diff + u_noise + u_adc + u_ds
            u_expand.append(u_stray_sys + self.k * u_1sigma)
            
        return [np.array(u_expand)[:, i] for i in range(np.array(u_expand).shape[1])]

