from setuptools import setup, find_packages

setup(
    name='tcc_python',
    version='1.0',
    description='A package for automating analyses with the TCC and the postprocessing of its output',
    author='Peter Crowther',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'pytest']
)
