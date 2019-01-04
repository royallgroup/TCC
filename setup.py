from setuptools import setup, find_packages

setup(
    name='tcc_python_scripts',
    version='1.0',
    description='A series of scripts for automating the TCC and postprocessing its output',
    author='Peter Crowther',
    packages=find_packages(),
    install_requires=['numpy', 'pandas']
)
