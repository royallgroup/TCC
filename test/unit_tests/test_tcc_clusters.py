"""Unit tests for detecting individual clusters."""

import pytest
import os
import numpy
import pandas
from glob import glob

from python_scripts.file_readers import xyz
from python_scripts.tcc import wrapper, structures

structures_to_test = glob('test/unit_tests/clusters/*.xyz')


def run_unit_test(cluster_path, bond_type):
    """Unit tests for clusters when running TCC.

    Failure if number of detected clusters does not match known composition of given structure.

    Args:
        cluster_path: path to isolated structure geometry to run TCC on
        bond_type: the type of bond to be used by the TCC
    """

    particle_coordinates = xyz.read(cluster_path).particle_coordinates
    cluster_name = os.path.split(cluster_path)[1].rstrip(".xyz")

    TCC = wrapper.TCCWrapper()
    if bond_type == "simple":
        TCC.input_parameters['Simulation']['bond_type'] = wrapper.BondType.simple
        TCC.input_parameters['Simulation']['rcutAA'] = 1.1

    tcc_result_report = TCC.run((100., 100., 100.), particle_coordinates)
    tcc_result_report = tcc_result_report['Number of clusters']

    report = pandas.DataFrame(tcc_result_report)
    report['Expected'] = 0
    for component, n in structures.composition[cluster_name].items():
        report['Expected'][component] = n

    try:
        assert numpy.all(report['Number of clusters'] == report['Expected'])
    except AssertionError:
        print(report)
        raise AssertionError from None


@pytest.mark.parametrize('path', structures_to_test)
def test_clusters(path):
    run_unit_test(path, "voronoi")
    run_unit_test(path, "simple")
