"""Unit tests for detecting individual clusters."""

import pytest
import sys
import os
import numpy, pandas
from glob import glob
from tcc import structures, xyz, wrapper

structures_to_test = glob('clusters/*.xyz')
sys.path += [os.path.abspath('../lib')]


def run_static_tcc_simple_bonds(x, box=(100., 100., 100.)):
    """Run TCC with pairwise cutoff bond detection algorithm.

    Args:
        x: cluster coordinates
        box: box side lengths for boundary conditions
    Returns:
        pandas table giving the static cluster information
    """
    TCC = wrapper.TCCWrapper()
    TCC.input_parameters['Simulation']['bond_type'] = wrapper.BondType.simple
    return TCC.run(box,x)


def run_static_tcc_voronoi(x, box=(100., 100., 100.)):
    """Run TCC with voronoi bond detection algorithm.

    Args:
        x: cluster coordinates
        box: box side lengths for boundary conditions
    Returns:
        pandas table giving the static cluster information
    """
    TCC = wrapper.TCCWrapper()
    return TCC.run(box,x)


def run_test(path, bond_type):
    """Unit tests for clusters when running TCC.

    Failure if number of detected clusters does not match known composition of given structure.

    Args:
        path: path to isolated structure geometry to run TCC on
        bond_type: the type of bond to be used by the TCC
    """
    particle_coordinates = xyz.read(path)
    cluster_name = path.split('/')[-1].split('.xyz')[0]

    if bond_type == "voronoi":
        tcc_result_report = run_static_tcc_voronoi(particle_coordinates)['Number of clusters']
    else:
        tcc_result_report = run_static_tcc_simple_bonds(particle_coordinates)['Number of clusters']
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
    run_test(path, "simple")
    run_test(path, "voronoi")
