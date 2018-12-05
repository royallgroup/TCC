# Docs README #

The documentation for the TCC is generated using the Sphinx Python package.

Documentation is written in reStructuredText format and compiled by the Sphinx package into a HTML source. At compilation time Sphinx also scrapes docstrings from the Python scripts and adds these to the documentation.

C docstrings are scraped by Doxygen to XML format and then imported into sphinx using the Breathe package (https://breathe.readthedocs.io/en/latest/) for Sphinx. You do not need to install Breathe, it is included in the extern folder in the TCC root directory.

Requrements to compile documentation
======================================

Install: Doxygen (http://www.doxygen.nl/)
Install: Sphinx (http://www.sphinx-doc.org/en/master/) - (conda install sphinx/pip install sphinx)


Compiling documentation
=========================

To recompile the documentation after a change, from the docs folder use the commands::

cd doxygen
doxygen
cd ..
sphinx-build -b html ./source ./html

Once a commit is made to the master branch the docs folder is mirrored to https://royallgroup.github.io/TCC/, providing an online version of the documentation.