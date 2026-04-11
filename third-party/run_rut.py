#!/bin/python3
# -*- coding: utf-8 -*-
"""
@author: aderu
"""

import os,sys
import json
import argparse
import shutil
from matplotlib import pyplot as plt
import numpy as np
import netCDF4 as nc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Source.S2Reader import S2Processor
from Source.S2RUT import S2RUT_L1

S2_BAND_NAMES = ["B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B10","B11","B12"]
S2_BAND_SAMPLING = {'B01':60, 'B02':10, 'B03':10, 'B04':10, 'B05':20, 'B06':20, 'B07':20, 'B08':10, 'B8A':20, 'B09':60, 'B10':60, 'B11':20, 'B12':20}

def json_none(param):
    if param=="None" or param=="":
        return None
    else:
        return float(param)
    

def run_s2rut_image(config):
    '''
    S2-L1C uncertainty computation for full image or AOI, without spectral error correlation.
    Work band per band to minimize RAM consumption.

    Parameters
    ----------
    config : str, path of config file

    Returns
    -------
    Save uncertainty product in reflectance dimension

    '''

    # check if inputs bands are not none
    if config["input_bands"] is None:
        bandlist = S2_BAND_NAMES
    else:
        bandlist =  config["input_bands"] 
        
    # check uncertainty type to output
    fname = os.path.basename(config['input_L1C']).split('.')[0]
    
    print("S2 RUT L1: image/AOI uncertainty")
    print("Processing:", os.path.basename(config['input_L1C']))
    print("Processing bands:", bandlist)
    do_plot = config['doplot'] 
    do_cont = config['unc_per_contributor']
    if do_cont:
        print("warning: computing unc contributors takes a long time!")

    # process individually each band to avoid RAM consumption in case of full image
    for i, band in enumerate(bandlist): 
        print(band)
        s2proc = S2Processor()
        s2proc.product = config['input_L1C'] 
        s2proc.noise_model = config['input_noise_model']
        s2proc.output = config['path_output']
        s2proc.selected_bands = [band]
        s2proc.lat_centre = json_none(config['roi_lat'])
        s2proc.lon_centre = json_none(config['roi_lon'])
        s2proc.w = json_none(config['roi_width'])
        s2proc.h = json_none(config['roi_height'])
        
        # check if roi size is larger than 1 pixel
        roi_pix_x = json_none(config['roi_width'])/S2_BAND_SAMPLING[band] 
        roi_pix_y = json_none(config['roi_height'])/S2_BAND_SAMPLING[band] 
        if roi_pix_x<2 or roi_pix_y<2:
            print("ERROR: ROI size is less than 2 pixels, skipping band")
            continue
            
        # retrieve L1C radiometry
        s2proc.get_data(get_va=False, doplot=False)
        img = s2proc.L1C_ref[0]

        # Compute uncertainty (absolute value, in reflectance dimension)
        RUTl1 = S2RUT_L1(config["input_contributors"])
        u_ref, u_cont = RUTl1.unc_calculation_abs(img, band, 0, s2proc.metadatadict, s2proc.sun_zenith, do_contributor=do_cont)
        
        # write output product
        newfile = os.path.join(s2proc.output, fname, band+"_unc_abs.jp2") 
        s2proc.write_unc(u_ref, newfile, band, s2proc.ds_profile[0])
        if do_cont:
            # write each contributor
            for u in u_cont.keys():
                print("  ", u)
                newfile = os.path.join(s2proc.output, fname, 'contributors', band+"_"+u+".jp2") 
                s2proc.write_unc(u_cont[u], newfile, band, s2proc.ds_profile[0])
    
        # copy input parameter for tracing configuration
        shutil.copy2(config["input_contributors"], os.path.dirname(newfile))
    
        # Plot
        if do_plot:
            fig1, ax1 = plt.subplots(1, 2, figsize=(8,5))
            plt.subplots_adjust(left=0.05, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.1)
            im0 = ax1[0].imshow(img, interpolation='None')
            ax1[0].set_title(band+' Reflectance')
            divider = make_axes_locatable(ax1[0])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cb = fig1.colorbar(im0, cax=cax, orientation='vertical', label='Reflectance')    
            cb.ax.tick_params(labelsize=8)
            
            img[img==0] = np.nan
            im1 = ax1[1].imshow(u_ref*100/img, interpolation='None', vmax=np.nanmean(u_ref*100/img)*3)
            ax1[1].set_title(band+' Uncertainties ')
            divider = make_axes_locatable(ax1[1])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cb = fig1.colorbar(im1, cax=cax, orientation='vertical', label="unc [%]")   
            cb.ax.tick_params(labelsize=8)

            for ax in ax1:
                ax.tick_params(axis='both', which='major', labelsize=8)
            fig1.savefig(os.path.join(s2proc.output, fname, band+"_unc_plot.png"))

        del u_ref, u_cont
        
    # Show all plots at once
    plt.show()
    print("Output files written in", config['path_output'])
    
    

def run_s2rut_spectralcorr(config):
    '''
    S2-L1C uncertainty computation with spectral error correlation.
    Uncetainties computed for all bands, but only for a single pixel.
    !!! NOT AVAILABLE YET !!!

    Parameters
    ----------
    config : str, path of config file

    Returns
    -------
    None.

    '''
       
    print("S2 RUT L1: spectral error correlation uncertainty")
    print("Only available for a single pixel : AOI size set to 1 pixel")
    do_plot = config['doplot'] 
    
    s2proc = S2Processor()
    s2proc.product = config['input_L1C'] 
    s2proc.noise_model = config['input_noise_model']
    s2proc.output = config['path_output']
    s2proc.lat_centre = json_none(config['roi_lat'])
    s2proc.lon_centre = json_none(config['roi_lon'])
    
    # hard coded parameters in this mode
    s2proc.selected_bands = ["B01","B02","B03","B04","B05","B06","B07","B08","B09","B10","B11","B12"]
    s2proc.w = 1
    s2proc.h = 1
    
    # Retrieve L1C radiometry
    s2proc.get_data(get_va=False, doplot=False)
       
    # Compute uncertainty
    RUTl1 = S2RUT_L1(config["input_contributors"])
    unc = []
    for i, band in enumerate(s2proc.selected_bands): 
        # u_ref, u_sig, u_sys, u_cont = RUTl1.unc_calculation(s2proc.L1C_ref[i], band, s2proc.metadatadict, s2proc.sun_zenith)
        u_ref = RUTl1.unc_calculation(s2proc.L1C_ref[i], band, i, s2proc.metadatadict, s2proc.sun_zenith)
        unc.append(u_ref)

    # Compute Spectral Error correlation
    nsamp = config['SEC_sample']
    u_corr = np.array(RUTl1.unc_spectralcorrelation(s2proc.L1C_rad, s2proc.metadatadict, nsamp))
    
    # write output
    fname = os.path.basename(config['input_L1C']).split('.')[0]
    new_file = os.path.join(config['path_output'], fname, 'spectral_err_corr_'+str(nsamp)+'.nc')
    s2proc.write_spectral_corr_error(unc, u_corr, new_file)

    if do_plot:
        plt.figure()
        plt.plot(np.squeeze(unc))
        plt.title("unc for given pixel")
        plt.xlabel('bands')
        plt.ylabel('unc [%]')
    
        plt.figure()
        for i in range(nsamp):
            plt.plot(u_corr[:,i])
        plt.title("unc spectral correlation per band")
        plt.xlabel('bands')
        plt.ylabel('unc [%]')

    # Show all plots at once
    plt.show()
    
    fname = os.path.join(config['path_output'], os.path.basename(config['input_L1C']).split('.')[0]+'.nc')
    f = nc.Dataset(fname, 'w', format='NETCDF4')
    f.createDimension('bands', len(s2proc.selected_bands) )
    f.createDimension('sample', config['SEC_sample'] )
    var1 = f.createVariable('pixel_unc', 'f4', ('bands') )
    var2 = f.createVariable('spectral_correlation_error', 'f4', ('bands', 'sample') )
    var1[:] = np.squeeze(unc)
    var2[:] = u_corr
    f.history ="S2-RUT-L1 spectrale correlation error"
    f.close()



#### MAIN
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                     description="S2-RUT-L1C tool")
    parser.add_argument("config_file",  help='path of RUT config file')
    args = parser.parse_args()


    # try:
    # Read config file
    with open(args.config_file) as f:
        config = json.load(f)
       
    # check if contributor file exists
    fcontrib = config["input_contributors"]
    if not os.path.exists(fcontrib):
        print('ERROR: uncertainty contributor file is missing')
        print('check Data/input_contributors.json in repository')
        sys.exit()   
       
    # launch RUT
    run_s2rut_image(config)
    
    # except Exception as error:
    #     print("An error occurred:", error) 
    #     exc_type, exc_obj, exc_tb = sys.exc_info()
    #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    #     print(exc_type, fname, exc_tb.tb_lineno)


