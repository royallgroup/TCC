# Docs README #

The documentation for the TCC is generated using the Sphinx Python package.

Documentation is written in reStructuredText format and compiled by the Sphinx package into a HTML source. At compilation time Sphinx also scrapes docstrings from the Python scripts and adds these to the documentation. Extraction of docstrings from the C source is planned but not yet in place.

To compile the documentation after a change use the command `sphinx-build -b html ./source ./html` from the docs folder.

Once a commit is made to the master branch the docs folder is mirrored to https://royallgroup.github.io/TCC/, thus providing an online version of the documentation.