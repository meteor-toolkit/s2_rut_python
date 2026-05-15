# s2_rut_python

Pure Python port of ESA's Sentinel-2 radiometric uncertainty tool.

`s2_rut_python` computes per-pixel radiometric uncertainties for Sentinel-2 Level 1C
reflectance data, following the methodology of Gorrono *et al.* (2017). Outputs are
stored in the `obsarray` uncertainty accessor format.

> **Warning:** This software is in beta. Results should be used with
> caution. Please share any feedback via the issue tracker.

## Usage

### Virtual environment

It is always recommended to use a virtual environment for each Python project.
Use your preferred environment manager, or create one with:

```bash
python -m venv venv
```

Activate it on Windows with `venv\Scripts\activate`, or on macOS/Linux with
`source venv/bin/activate`.

### Installation

Install the package and its core dependencies:

```bash
pip install -e .
```

Optional extras are available depending on your use case:

```bash
pip install -e ".[dev]"      # Development tools (ruff, mypy, pytest, …)
pip install -e ".[docs]"     # Documentation build (sphinx, …)
```

### Development

Install the pre-commit hooks after cloning:

```bash
pre-commit install
```

When you commit, `ruff` will lint and format your code. If it makes
corrections the commit will be aborted so you can review the changes — just
commit again once you are happy.

Run the test suite with:

```bash
pytest
```

## Compatibility

`s2_rut_python` requires Python 3.11 or later and is tested on Python 3.11, 3.12,
and 3.13.

## Licence

`s2_rut_python` is released under the GNU Lesser General Public License v3 (LGPLv3).
See the [LICENSE](LICENSE) file for the full licence text.

## Authors

`s2_rut_python` is developed and maintained by the
[MetEOR Toolkit Team](mailto:team@comet-toolkit.org).

## Citations

Gorrono, J.; Fomferra, N.; Peters, M.; Gascon, F.; Underwood, C.I.; Fox, N.P.;
Kirches, G.; Brockmann, C. A Radiometric Uncertainty Tool for the Sentinel 2 Mission.
*Remote Sens.* 2017, 9, 178. <https://doi.org/10.3390/rs9020178>
