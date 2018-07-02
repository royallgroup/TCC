"""Unit tests for detecting individual clusters."""

import pytest
import os
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

    particle_coordinates = list(xyz.read(cluster_path, 1))[0].particle_coordinates
    cluster_name = os.path.split(cluster_path)[1].rstrip(".xyz")

    tcc_parameters = wrapper.TCCWrapper()
    if bond_type == "voronoi_short":
        tcc_parameters.input_parameters['Simulation']['voronoi_parameter'] = 0.82
    elif bond_type == "voronoi_long":
        tcc_parameters.input_parameters['Simulation']['voronoi_parameter'] = 1
    else:
        print("Unknown bond type.")
        raise TypeError()

    tcc_result_report = tcc_parameters.run((100., 100., 100.), particle_coordinates)

    report = pandas.DataFrame(tcc_result_report['Number of clusters'])

    try:
        if bond_type == "voronoi_short":
            assert report['Number of clusters'][cluster_name] == structures.voronoi_short_clusters[cluster_name]
        elif bond_type == "voronoi_long":
            assert report['Number of clusters'][cluster_name] == structures.voronoi_long_clusters[cluster_name]
        else:
            print("Unknown bond type.")
            raise TypeError()
    except AssertionError:
        print(report)
        raise AssertionError from None


@pytest.mark.parametrize('path', structures_to_test)
def test_voronoi_long_clusters(path):
    run_unit_test(path, "voronoi_long")


@pytest.mark.parametrize('path', structures_to_test)
def test_voronoi_short_clusters(path):
    run_unit_test(path, "voronoi_short")
