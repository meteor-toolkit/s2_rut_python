import s2_rut

def run_RUT(product,params):
    """
    product is an eopy instance
    params  is a dictionary with the set of parameters to run the RUT
            including:
                all_flags                                  Sets all the flags to True
                all_contribs                               Runs the RUT with default settings - if none of the flags (this or those below) are set
                                                           This flag is set to True.
                PADC_quantisation=<boolean>                Analog to Digital conversion at the Video Chain Unit on-board the MSI
                                                           Default value is 'False'.
                Pband_names=<string,string,string,...>     The bands for which the uncertainty shall be computed
                                                           Default value is 'B1,B2,B3,B4,B5,B6,B7,B8,B8A,B9,B10,B11,B12'.
                Pcoverage_factor=<double>                  The value of the coverage factor for the uncertainty evaluation (
                                                           Default value is '1.0'.
                PCrosstalk=<boolean>                       Focal plane (Optical) and Front-End Electronics (electrical) inter-band signal
                                                           Default value is 'False'.
                PDiffuser-absolute_knowledge=<boolean>     Knowledge on the diffuser reflectance factor (BRF)
                                                           Default value is 'False'.
                PDiffuser-cosine_effect=<boolean>          Cosine correction knowledge as a consequence of angular noise
                                                           Default value is 'False'.
                PDiffuser-straylight_residual=<boolean>    Residual for the correction of the stray-light during in-flight diffuser calibration
                                                           Default value is 'False'.
                PDiffuser-temporal_knowledge=<boolean>     Estimated effect of the diffuser degradation in space
                                                           Default value is 'False'.
                PDS_stability=<boolean>                    Residual thermal fluctuations of the detector offset along the orbit
                                                           Default value is 'False'.
                PGamma_knowledge=<boolean>                 Knowledge on te correction for non-linearity and non-uniformity
                                                           Default value is 'False'.
                PInstrument_noise=<boolean>                Noise (shot, thermal, etc.) introduced by the silicon and CMT detectors of the MSI instrument
                                                           Default value is 'False'.
                PL1C_image_quantisation=<boolean>          Effect of the finite resolution of the L1C reflectance factor
                                                           Default value is 'False'.
                POOF_straylight-random=<boolean>           Focal plane Out-Of-Field light that results in a random spatial dispersion
                                                           Default value is 'False'.
                POOF_straylight-systematic=<boolean>       Telescope Out-Of-Field light that results in a positive bias
                                                           Default value is 'False'.
    """
    # parse the commands input in the dictionary
    param_keys = params.keys()
    gpt_args   = []

    # all contributions flag
    if "all_contribs" in param_keys:
        if params["all_contribs"]:
            gpt_args = []

    else:
        # more effective way to do this but whatever
        args = ["PADC_quantisation",
                "PCrosstalk",
                "PDiffuser-absolute_knowledge",
                "PDiffuser-cosine_effect",
                "PDiffuser-straylight_residual",
                "PDiffuser-temporal_knowledge",
                "PDS_stability",
                "PGamma_knowledge",
                "PInstrument_noise",
                "PL1C_image_quantisation",
                "POOF_straylight-random",
                "POOF_straylight-systematic"]
        gpt_args = populate_args(params,[],args)
        if "Pband_names" in param_keys:
            gpt_args.append("-Pband_names="+params["Pband_names"])
        else:
            pass
        if "Pcoverage_factor" in param_keys:
            gpt_args.append("-Pcoverage_factors="+params["Pcoverage_factor"])
        else:
            pass
        if "all_flags" in param_keys:
            for arg in args:
                if ("-"+arg+"=False" in gpt_args) or ("-"+arg+"=True" in gpt_args):
                    gpt_args.remove("-"+arg+"=False")

    # create temporary directory to store RUT files in and move there
    os.mkdir(mdir)
    cdir = os.getcwd()
    os.chdir(mdir)

    ## save current product to a file located in mdir with filename saved in "temp_fname"
    # Sam: you should know the best way to do this...
    product.save(mdir+"preProd.dim") # must save as a snappy product

    # run the RUT from the command line
    sp.call([gpt,"gpt S2RutOp ",mdir+temp_fname+" "," ".join(gpt_args)],shell=False)

    # read the RUT data in
    # standard eopy read in should work here - Sam: this doesn't work - it doesn't appear that it likes .dim files
    # this should just go down the standard MTD...xml route?
    openRUT = eopy.Product("target.dim")
    ## attribute array stuff here that you need to do
    ## not sure what exactly you want to do with the info

    # close the product
    openRUT.closeProduct()

    # clean up temporary directories and files
    os.chdir(cdir)
    os.remove(mdir+"target.dim")
    os.remove(mdir+"preProd.dim")
    shutil.rmtree(mdir+"target.data")


def populate_args(param_dict,arg_list,args):
    for arg in args:
        if arg in param_dict.keys():
            if param_dict[arg]:
                args_list.append("-%s=True"%arg)
            else:
                args_list.append("-%s=False"%arg)
        else:
            args_list.append("-%s=False"%arg)

    return args_list
