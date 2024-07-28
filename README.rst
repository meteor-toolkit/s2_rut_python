s2_rut_python
=============

Pure Python port of Javi's Sentinel-2 radiometric uncertainty tool. Utilises the s2_rut_algo.py file containing the
individual uncertainty parameters and the calculations for finding the individual uncertainty contributions. Also,
s2_l1_rad_conf.py file is required, which contains band- and spacecraft-specifc values required for calculating the
uncertainties.

Usage
=====
Files are read into an xarray format for use with our other tools. Eoio capabilities are used for accessing the metadata.

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

`s2_rut_python` was written by `Sam Hunt <sam.hunt@npl.co.uk>`_.
