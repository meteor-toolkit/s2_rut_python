s2_rut_python
=============

Pure Python portion of Javi's Sentinel-2 radiometric uncertainty tool.

Usage
=====
Files are read into an xarray format for use with our other tools. eoio capabilities are used for accessing and
extracting the metadata, used for calculating the uncertainties.

Band selection takes one band or a list of **Sentinel-2 band names**: B01, B02, B03, B04, B05, B06, B07, B08, B8A, B09,
B10, B11, B12.

The current interface entry point is ``S2RUTTool`` and uncertainties can be returned either:

* grouped by correlation component (``u_systematic_<band>``, ``u_random_<band>``), or
* per contributor (for example ``u_noise_<band>``).

Each individual type of components is categorised according to Gorrono *et al.* (2017).


Usage of this module is as follows ::

    from s2_rut_python.interface import S2RUTTool

    s2rut = S2RUTTool()
    ds_out = s2rut.run(
        ds=s2_ds,
        data_vars=['B01'],
        group_unc=True,
    )

To run per-contributor uncertainties ::

    ds_out = s2rut.run(
        ds=s2_ds,
        data_vars=['B01'],
        group_unc=False,
    )

To limit contributors to a subset ::

    ds_out = s2rut.run(
        ds=s2_ds,
        data_vars=['B01'],
        group_unc=False,
        subset_unc=['noise', 'adc', 'gamma'],
    )

Uncertainty variables are stored in the obsarray uncertainty accessor (for example
``ds_out.unc['B01']``). Reflectance zero-values are masked to ``NaN`` before return.

Example
=======
An up-to-date runnable example is available in ``examples/run_example.py``.
It demonstrates reading Sentinel-2 L1C data with eoio, selecting an ROI subset for the current product, running ``S2RUTTool.run(...)``,
and plotting reflectance, random uncertainty [%], and systematic uncertainty [%].

Virtual environment
-------------------

It's always recommended to make a virtual environment for each of your python
projects. Use your preferred virtual environment manager if you want and
activate it for the rest of these commands. If you're unfamiliar, read
https://realpython.com/python-virtual-environments-a-primer/. You can set one up
using::

    python -m venv venv

and then activate it on Windows by using ``venv/Scripts/activate``. 

Installation
------------

Install your package and its dependencies by using::

    pip install -e .

Development
-----------

For developing the package, you'll want to install the pre-commit hooks as well. Type::

    pre-commit install


Note that from now on when you commit, `black` will check your code for styling
errors. If it finds any it will correct them, but the commit will be aborted.
This is so that you can check its work before you continue. If you're happy,
just commit again. 

Compatibility
-------------

Licence
-------

Authors
-------

s2_rut_python has been developed by `Rasma Ormane <rasma.ormane@npl.co.uk>`_. `Sam Hunt <sam.hunt@npl.co.uk>`_. and `Maddie Stedman <maddie.stedman@npl.co.uk>`_.

Citations
-------

Gorrono, J.; Fomferra, N.; Peters, M.; Gascon, F.; Underwood, C.I.; Fox, N.P.; Kirches, G.; Brockmann, C.
A Radiometric Uncertainty Tool for the Sentinel 2 Mission. Remote Sens. 2017, 9, 178. https://doi.org/10.3390/rs9020178