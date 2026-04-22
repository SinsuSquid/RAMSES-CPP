# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'RAMSES-CPP'
copyright = '2026, RAMSES-CPP Contributors'
author = 'RAMSES-CPP Contributors'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_simplepdf',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# SimplePDF settings (used by .readthedocs.yaml)
simplepdf_vars = {
    'primary': '#333333',
    'primary_opaque': '#33333355',
    'secondary': '#0066cc',
    'cover': '#ffffff',
    'white': '#ffffff',
    'links': '#0066cc',
}
