#!/usr/bin/env python
#
# s2_rut_python documentation build configuration file
#

from importlib.metadata import version as _pkg_version, PackageNotFoundError

# Private NPL packages not available in the plain docs-build environment.
# autodoc_mock_imports makes these importable as stubs during autodoc processing.
_MOCK_MODULES = [
    "obsarray",
    "processor_tools",
    "processor_tools.config",
    "processor_tools.context",
]

try:
    _version = _pkg_version("s2_rut_python")
except PackageNotFoundError:
    _version = "0.0.0"

project_title = "s2_rut_python".replace("_", " ").title()


# -- General configuration ---------------------------------------------

default_role = "code"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_design",
]

templates_path = ["_templates"]

source_suffix = ".rst"

master_doc = "index"

project = project_title
copyright = "MetEOR Toolkit Team"
author = "MetEOR Toolkit Team"

version = _version
release = _version

language = "en"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_mock_imports = _MOCK_MODULES

pygments_style = "sphinx"

todo_include_todos = False


# -- Options for HTML output -------------------------------------------

html_theme = "sphinx_book_theme"

html_static_path = ["_static"]

htmlhelp_basename = "s2_rut_pythondoc"


# -- Options for LaTeX output ------------------------------------------

latex_elements = {}

latex_documents = [
    (
        "content/user/user_guide",
        "user_manual.tex",
        "{}: User Guide".format(project_title),
        "MetEOR Toolkit Team",
        "manual",
    ),
    (
        "content/user/atbd",
        "atbd.tex",
        "{}: Algorithm Theoretical Basis Document".format(project_title),
        "MetEOR Toolkit Team",
        "manual",
    ),
]
