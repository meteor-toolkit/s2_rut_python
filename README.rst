s2_rut_python
=============

Pure Python portion of Javi's Sentinel-2 radiometric uncertainty tool.

Usage
=====
Files are read into an xarray format for use with our other tools. Eoio capabilities are used for accessing  and
extracting the metadata, used for calculating the uncertainties.

User specifies three components::

    s2rut=S2RUT()
    test = s2rut.run_rut(
        data_set=xarray_dataset
        band_names=['list of band names']
        unc_info='output type'
       )

Band names takes in one band or a list of **Sentinel-2 band names**: B01, B02, B03, B03, B04, B05, B06, B07, B08, B8A, B09,
B10, B11, B12.

There are two types of possible uncertainty information types: **'total'** which returns the total uncertainty, while
**'components'** provides three separate uncertainties: **'systematic**, **'random'**, and **'structured'**.

Each individual type of components is categorised according to Gorrono *et.al.* (2017).




Example
=======
An example illustrating the circumstances in which your code would normally be used is available in run_s2rut_example.py.
Two different input cases are demonstrated: Sentinel-2A and Sentinel-2B in .SAFE (Standard Archive Format for Europe)
data format.




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

Install your package and its dependancies by using::

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

s2_rut_interface was written by `Rasma Ormane <rasma.ormane@npl.co.uk>`_.

s2_rut_python` was written by `Sam Hunt <sam.hunt@npl.co.uk>`_.

Citations
-------

Gorrono, J.; Fomferra, N.; Peters, M.; Gascon, F.; Underwood, C.I.; Fox, N.P.; Kirches, G.; Brockmann, C.
A Radiometric Uncertainty Tool for the Sentinel 2 Mission. Remote Sens. 2017, 9, 178. https://doi.org/10.3390/rs9020178