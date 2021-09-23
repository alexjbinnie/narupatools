# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../narupatools/src"))

# -- Project information -----------------------------------------------------

project = "narupatools"
copyright = "2021, Alex Jamieson-Binnie"
author = "Alex Jamieson-Binnie"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.doctest",
    "nbsphinx",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"
html_theme_options = {
    "font_family": "Lato",
    "font_size": "16px",
    "head_font_family": "Roboto Slab",
    "page_width": "1000px",
    "fixed_sidebar": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = [
    "style.css",
]

autodoc_mock_imports = ["lammps"]

autosummary_generate = True

nbsphinx_execute = "always"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase", None),
    "mdanalysis": ("https://docs.mdanalysis.org/1.0.0", None),
    "mdtraj": ("https://mdtraj.org/1.9.4", None),
    "narupa": ("https://narupa.readthedocs.io/en/latest", None),
    "openmm": ("http://docs.openmm.org/7.1.0/api-python/", None),
}

set_type_checking_flag = True

inheritance_graph_attrs = dict(rankdir="TB")
