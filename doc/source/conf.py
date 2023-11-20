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
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'chemcomp'
copyright = '2021, Aaron Schneider'
author = 'Aaron Schneider'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax', 'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.extlinks',
    'sphinx.ext.viewcode',
    'numpydoc',
    'nbsphinx',
    'recommonmark',
    'sphinx.ext.autosectionlabel'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# https://stackoverflow.com/questions/12206334/sphinx-autosummary-toctree-contains-reference-to-nonexisting-document-warnings
# numpydoc_show_class_members = False

# NOT to sort autodoc functions in alphabetical order
autodoc_member_order = 'bysource'
autosummary_generate = True
numpydoc_show_class_members = False

# To avoid installing xESMF and all its dependencies when building doc
# https://stackoverflow.com/a/15912502/8729698
autodoc_mock_imports = ['xarray', 'scipy', 'ESMF', 'xesmf', 'xmitgcm', 'ESMPy', 'xgcm', 'matplotlib']

# avoid automatic execution for notebooks
nbsphinx_execute = 'never'
# Numpydoc settings



# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
