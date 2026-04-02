.. _quickstart:

################
Quickstart Guide
################

Sentinel-2 Radiometric Uncertainty Tool
+++++++++++++++++++++++++

``s2_rut_python`` provides a Python interface around ESA's Sentinel-2 Radiometric Uncertainty
Tool (S2-RUT) to compute pixel-wise radiometric uncertainty from Sentinel-2 L1C reflectance data.

The main entry point is :py:class:`s2_rut_python.interface.S2RUTTool`.

Theoretical background and original SNAP implementation references are listed in :ref:`Citations <citations>` below.

Requirements
------------

Input must be an :py:class:`xarray.Dataset` containing:

* Sentinel-2 reflectance band(s) (``B01`` ... ``B12``)
* solar zenith angles in either:

    * ``solar_zenith_angle`` with matching band shape, or
    * ``solar_zenith_angle_interp`` with matching band shape
* dataset metadata:

    * ``platform``
    * ``quantification_level``
    * ``reflectance_conversion_u``
* per-band metadata in ``ds[band].attrs['product_metadata']``:

    * ``solar_irradiance``
    * ``noise_model_alpha``
    * ``noise_model_beta``
    * ``physical_gains``
    * ``radiometric_offset``

Basic Usage
-----------

Run grouped uncertainties (default) for selected bands:

::

        from s2_rut_python.interface import S2RUTTool

        s2rut = S2RUTTool()
        ds_out = s2rut.run(
                ds=ds,
                data_vars=["B01", "B09"],
                group_unc=True,
        )

Run per-contributor uncertainties:

::

        ds_out = s2rut.run(
                ds=ds,
                data_vars=["B01"],
                group_unc=False,
        )

Limit contributors to a subset:

::

        ds_out = s2rut.run(
                ds=ds,
                data_vars=["B01"],
                group_unc=False,
                subset_unc=["noise", "adc", "gamma"],
        )

Run Parameters
--------------

``S2RUTTool.run(ds, data_vars=True, group_unc=True, subset_unc=None)``

* ``data_vars``:

    * ``True``: process all Sentinel-2 bands listed in ``MEAS_VAR_RES``
    * ``str``: process one band (for example ``"B01"``)
    * ``list[str]``: process explicit band list

* ``group_unc``:

    * ``True``: returns grouped uncertainties by correlation component

        * ``u_systematic_<band>``
        * ``u_random_<band>``

    * ``False``: returns per-contributor uncertainties (for example ``u_noise_<band>``)

* ``subset_unc``:

    * ``None``: include all contributors
    * list of contributor names from:

        * ``noise``
        * ``stray_sys``
        * ``stray_rand``
        * ``xtalk``
        * ``adc``
        * ``ds``
        * ``gamma``
        * ``diff_abs``
        * ``diff_temp``
        * ``diff_cos``
        * ``diff_sl``
        * ``ref_quant``
        * ``geoloc``

Output
------

The output is the input dataset with uncertainty variables stored in the obsarray uncertainty accessor
(``ds_out.unc[band]``) for each processed band.

Each uncertainty variable includes metadata such as:

* ``long_name``
* ``description``
* ``units`` (matching the source reflectance band)
* ``standard_name``
* ``pdf_shape``
* error-correlation metadata via ``err_corr``

Zero Reflectance Handling
-------------------------

Pixels where reflectance equals zero are masked to ``NaN`` in the reflectance band before return.
Uncertainty arrays are already masked consistently during creation.

Example Script
--------------

An up-to-date runnable example is provided in:

::

        examples/run_s2rut_example.py

It demonstrates:

* reading Sentinel-2 L1C data with ``eoio``
* running ``S2RUTTool.run``
* plotting reflectance, random uncertainty [%], and systematic uncertainty [%]

.. _citations:

Citations
---------

Gorrono, J.; Fomferra, N.; Peters, M.; Gascon, F.; Underwood, C.I.; Fox, N.P.; Kirches, G.; Brockmann, C.
A Radiometric Uncertainty Tool for the Sentinel 2 Mission. Remote Sens. 2017, 9, 178. https://doi.org/10.3390/rs9020178

S2-RUT. STEP. ESA. 2017. https://step.esa.int/main/wp-content/help/versions/9.0.0/snap-community-plugins/org.esa.snap.snap.rut/S2RutDocumentation.html



