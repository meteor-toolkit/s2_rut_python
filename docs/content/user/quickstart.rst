.. _quickstart:

################
Quickstart Guide
################

Sentinel-2 Radiometric Uncertainty Tool
+++++++++++++++++++++++++

S2-RUT is a Radiometric Uncertainty tool for the Sentinel 2 mission. Theoretical breakdown
of the tool, specifying individual uncertainty contributions, algorithms and equations used for calculations, and
original documentation is available at Gorrono `et.al.` (2017) and S2-RUT (2017).

It is available as a `SNAP` (Sentinel Application Platform) plugin, for obtaining pixel-by-pixel uncertainties for
Sentinel 2A and 2B products.

On SNAP, once a S2 product is selected, user can define the bands of interest (B01-B12), the desired coverage factor
(1.0, 2.0, or 3.0), and uncertainty contributions of interest:

* Instrument noise
* Out-of-Field straylight-systematic
* Out-of-Field straylight-random
* Crosstalk
* ADC quantisation
* DS stability
* Gamma knowledge
* Diffuser-absolute knowledge
* Diffuser-temporal knowledge
* Diffuser-cosine effect
* Diffuser-straylight residual
* L1C image quantisation

The plugin then uses the **s2_rut_algo** module to calculate each of these contributions, based on the uncertainty parameters
retrieved from the product and its metadata. There are 12 parameters:

* a - physical gains
* e_sun -
* u_sun - solar irradiance
* tecta - solar zenith angle (interpolated)
* quant - efect of the finite resolution of the L1C reflectance factor
* alpha - noise
* beta - noise
* u_diff_temp - estimated effect on the diffuser degeneration in space
* u_diff_cos - cosine correction knowledge as a consequence of angular noise
* u_diff_k - residual for the correction of the stray-light during in-flight diffuser calibration
* u_ADC - analogue to digital conversion at the Video Chain Unit on-board the MSI
* u_gamma - knowledge on the correction for non-linearity and non-uniformity

Depending on the specified uncertainty contributions, the relevant parameters are used to calculate the total radiometric
uncertainty which is returned in the form of a 2-dimensional array, of the same shape as the chosen band(s):

* B01, B09, B10 (60m res): (1830, 1830)
* B02, B03, B04, B08 (10m res): (10980, 10980)
* B05, B06, B07, BA8, B11, B12 (20m res): (5490, 5490)

Github containing source code for S2-RUT (written in python and Java) is availabe here: https://github.com/senbox-org/snap-rut/tree/master.

s2_rut_interface and SNAP overlaps
+++++++++++++++++++++++++

This purely python version of S2-RUT utilises only two modules: **s2_rut_algo** (uncertainty calculation algorithm) and **s2_l1_rad_conf**
(spacecraft-specific parameter values constant across all products, but varying for individual bands). The following values
are directly obtained from each of the modules:

**s2_rut_algo:**

* u_diff_cos: 0.4
* u_diff_k: 0.3
* u_ADC: 0.5
* u_gamma: 0.4

Note: these values were last checked on 01.08.2024.

**s2_l1_rad_conf:**

* u_diff_temp_rate - the rate of estimated effect on the diffuser degeneration in space


If the magnitude of either of these constant parameters changes in the source code, **s2_rut_interface** module
will adjust accordingly without having to change them manually.


s2_rut_interface vs. SNAP plugin
+++++++++++++++++++++++++

Uncertainty estimates of the  **s2_rut_interface** module have been successfully validated with SNAP outputs for multiple
S2 products. It has been observed that >97% of the product uncertainty outputs are identical for both S2-RUT application
methods. The remaining few pixels, however, appear to have dissimilar uncertainties, always having a difference of +/- 1.

Discrepancies
==============

Small deviations were noticed for randomly distributed pixels across the uncertainty outputs, when reading in products
eoio and performing uncertainty analysis using **s2_rut_interface** versus reading in the data on SNAP and utilising the
S2-RUT plugin tool.

Initially, the reason for was thought to be due to solar zenith angle interpolation methods used: **linear** for eoio
and **nearest** for SNAP.

The **s2_rut_interface** module provides the `S2RUT` class for calculating radiometric uncertainty in Sentinel-2 data.
This class interfaces with the `MyS2RUTAlgo` class, which extends functionality from `S2RutAlgo` in the `s2_rut_algo`
module.

Utilises the s2_rut_algo.py file containing the individual uncertainty parameters and the calculations for finding the individual uncertainty contributions. Also,
s2_l1_rad_conf.py file is required, which contains band- and spacecraft-specifc values required for calculating the uncertainties.

Input Parameters
-----------------

There are three input parameters that are required to run the s2_rut_interface module:

* ds: satellite data from a satellite product (Sentinel-2A or Sentinel-2B)
* band_names: desired bands

ds
^^^^^^^^^^^^^

Satellite data set must be read into an :py:class:`xarray.Dataset`. All bands of interest must be read in, before calling the **s2_rut_interface** module.
For each band, relevant metadata and parameter values are extracted and compiled in a dictionary in the **get_band_unc_parameters** method.

`Note: currently, module can accurately extract metadata from data that is read in using eoio. Xarrays extracted from SNAP,
of the same satellite data products have different names/chins for accessing metadata.`


band_names
^^^^^^^^^^^^^

Band names are provided as a list, where either one band can be analysed or all 13 together. The options are::

    band_names = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']


Constants
-----------------

At the beginning of the module, multiple lists of constants are defined, which are not xpected to change over time.

* **MEAS_VAR_RES** - list out all the bands and their respective resolutions (10, 20, or 60 m)
* **TIME_INIT** - specifies the date and time of the Sentinel-2A and -2B sensor launches. Used to calculate **u_diff_temp**.
* **CORR_TYPES** - categorises correlations (systematic, random, structured) for all possible uncertainty contributions
* **U_CONTRIBUTIONS** - list of all possible uncertainty contributions


Classes
-----------------

This module consists of two classes:

* **S2RUT** - the main class, with three methods (run_rut, return_unc_correlations, and get_band_unc_parameters)
* **MyS2RUTAlgo** - copy of a class (S2RutAlgo) from the original code (**s2_rut_algo** module). Overwrites the defined uncertainty
parameter values so they can be replaced with those extracted using eoio. Overwrites the **unc_select** instance, to enable
selecting individual uncertainty contributions.

Methods
-----------------

run_rut
^^^^^^^^

Runs the Sentinel 2 radiometric uncertainty tool, for selected bands and uncertainty types. Returns an xarray dataset with uncertainty variables for each specified band, containing information about its dimensions, coordinates and data variables::

    <xarray.Dataset>
    Dimensions:                    (y_60m: 1830, x_60m: 1830, x_5000m: 23,
                                    y_5000m: 23)
    Coordinates:
      * x_60m                      (x_60m) float64 15kB 7e+05 7e+05 ... 8.097e+05
      * y_60m                      (y_60m) float64 15kB 3.2e+06 3.2e+06 ... 3.09e+06
        ...
    Data variables:
        B01                        (y_60m, x_60m) float32 13MB 0.2499 ... 0.2225
        ...                         ...
        u_systematic_B01           (y_60m, x_60m) uint8 3MB 22 22 22 22 ... 22 22 22
        u_random_B01               (y_60m, x_60m) uint8 3MB 5 5 5 5 5 ... 5 5 5 5 5

return_unc_correlations
^^^^^^^^^^^^^^^^^^^^^^^

Returns dictionary of uncertainty correlations: systematic, random, and structured. For each type, the full contributions list
is available, where respective uncertainty parameters are classified as True or False, depending on their correlation status.

* **systematic**: Out-of-Field straylight-systematic, Crosstalk, Diffuser-straylight residual, Diffuser-temporal knowledge, Diffuser-absolute knowledge, DS stability, Diffuser-cosine_effect
* **random**: Out-of-Field straylight-random, Instrument noise,  ADC quantisation, Gamma knowledge, L1C image quantisation


This characterisation is in accordance to Gorrono `et.al.` (2017).

`Note: in the original **s2_rut_algo**, all uncertainty contributions except for **u_diff_temp** are used when calculating the total uncertainty. Hence if all individual uncertainties (random, systematic and structured) get added together in quadrature, the 'total' output versus a combination of all three slightly differs.`


get_band_unc_parameters
^^^^^^^^^^^^^^^^^^^^^^^

Extracts band-specific uncertainty parameters from the provided eoio data_set.

Example
-------
An example of the module used to calculate the uncertainty components of the following Sentinel-2A product::

    T:\ECO\EOServer\data\satellite\S2A_MSI\L1\PICS\Libya4\S2A_MSIL1C_20190103T090351_N0207_R007_T34RGS_20190103T110458.SAFE


is visible here::

    s2rut = S2RUT()
    test = s2rut.run_rut(
        ds=s2_l1c_ds,
        band_names=['B01', 'B09'],
    )

Running this code, gives the following output ::

    <xarray.Dataset>
    Dimensions:                    (y_60m: 1830, x_60m: 1830, x_5000m: 23,
                                    y_5000m: 23)
    Coordinates:
      * x_60m                      (x_60m) float64 15kB 7e+05 7e+05 ... 8.097e+05
      * y_60m                      (y_60m) float64 15kB 3.2e+06 3.2e+06 ... 3.09e+06
        ...
    Data variables:
        B01                        (y_60m, x_60m) float32 13MB 0.2499 ... 0.2225
        B09                        (y_60m, x_60m) float32 13MB 0.1716 ... 0.1651
        ...                         ...
        u_systematic_B01           (y_60m, x_60m) uint8 3MB 22 22 22 22 ... 22 22 22
        u_random_B01               (y_60m, x_60m) uint8 3MB 5 5 5 5 5 ... 5 5 5 5 5
        u_systematic_B09           (y_60m, x_60m) uint8 3MB 12 12 12 12 ... 12 12 12
        u_random_B09               (y_60m, x_60m) uint8 3MB 23 23 23 24 ... 24 23 23


Citations
-------

Gorrono, J.; Fomferra, N.; Peters, M.; Gascon, F.; Underwood, C.I.; Fox, N.P.; Kirches, G.; Brockmann, C.
A Radiometric Uncertainty Tool for the Sentinel 2 Mission. Remote Sens. 2017, 9, 178. https://doi.org/10.3390/rs9020178

S2-RUT. STEP. ESA. 2017. https://step.esa.int/main/wp-content/help/versions/9.0.0/snap-community-plugins/org.esa.snap.snap.rut/S2RutDocumentation.html



