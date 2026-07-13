.. _developer_guide:

################
Developer Guide
################

This guide describes how to set up a development environment and contribute to
``s2_rut_python``.

Installation
------------

Clone the repository and install the package with its development dependencies:

::

        pip install -e ".[dev]"

Install the pre-commit hooks after cloning:

::

        pre-commit install

When you commit, ``ruff`` will lint and format your code. If it makes corrections
the commit will be aborted so you can review the changes — just commit again once
you are happy.

Tests
-----

Run the test suite with:

::

        pytest

Static analysis
----------------

Run ``ruff`` and ``mypy`` before submitting changes:

::

        ruff format s2_rut_python
        ruff check s2_rut_python
        mypy s2_rut_python

Documentation
-------------

Build the documentation locally with:

::

        sphinx-build -b html docs docs/_build/html
