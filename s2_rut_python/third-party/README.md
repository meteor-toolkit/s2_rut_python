# S2-RUT-L1C
Stand alone tool computing Sentinel2-MSI L1C uncertainties. Developed by OPT-MPC, based on the <a href='https://github.com/senbox-org/snap-rut'>S2-RUT SNAP plugin</a>


## Installation & requirement
Python version 3.8 (or above) is required and the necessary python packages are listed in the environment.yml file. 
With <a href='https://docs.anaconda.com/miniconda/miniconda-install/'>conda</a> a dedicated S2RUT python environment can be setup with:
```
conda env create -f environment.yml
```

## Usage and configuration
The S2RUT run options are given in a json configuration file passed as an argument of the python fonction.
```
python3 run_rut.py config.json
```
The config.json file arguments are organised as follow:

**inputs parameters (required)**
- input_L1C: path to S2 L1C product (SAFE folder)
- input_bands: Choice of band to process, ex: ["B01", "B09"].
- input_contributors: path to uncertainty input contributors file.
- input_noise_model: path to noise model file, if "None" metadata will be used.

**ouputs parameters (required)**
- path_output: path of the output directory.
- unc_per_contributor: true or false, activate additional ouput file per contribuor. (Takes long time!!!)
- doplot: true or false, activate the plotting of the AOI. Plots are saved in the output folder.

**ROI definition (optional)**
- roi_lat: Centre latitude (decimal) of the ROI, if "None" the centre of image will be used.
- roi_lon: Centre longitude (decimal) of the ROI, if "None" the centre of image will be used. 
- roi_width: Width in meters of the ROI, If "None" the full image is used.
- roi_height: Height in meters of the ROI, If "None" the full image is used.


**output files**  
Absolute per pixel uncertainty in reflectance dimension are saved in jpeg2000 file, similarly to L1C radiometry (one file per band).


## Authors:
Alexis Deru @ ACRI-ST  
Sébastien Clerc @ ACRI-ST  
Javier Gorroño @ Universitat Politècnica de València  

