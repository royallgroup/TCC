"""Unit tests for detecting individual clusters."""

import pytest

import sys, os
sys.path.append(os.path.abspath('../lib'))
from tcc import structures, xyz, wrapper

def run_static_tcc_simple_bonds(x, box=[100.,100.,100.]):
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

def run_static_tcc_voronoi(x, box=[100.,100.,100.]):
    """Run TCC with voronoi bond detection algorithm.

    Args:
        x: cluster coordinates
        box: box side lengths for boundary conditions
    Returns:
        pandas table giving the static cluster information
    """
    TCC = wrapper.TCCWrapper()
    return TCC.run(box,x)

@pytest.mark.parametrize('cluster', structures.clusters)
def test_simple_bonds(cluster):
    """Unit tests for clusters when running TCC with simple pairwise cutoff bond detection algorithm.

    Failure if number of detected clusters does not match known composition of given structure.

    Args:
        cluster: isolated structure to run TCC
    """
    try: x = xyz.read('clusters/%s.xyz' % cluster)
    except: return # temporary until we have sample coordinates for all known structures!

    found = run_static_tcc_simple_bonds(x)['Number']
    for component,n in structures.composition[cluster].items():
        assert n == int(found[component])

@pytest.mark.parametrize('cluster', structures.clusters)
def test_voronoi(cluster):
    """Unit tests for clusters when running TCC with simple pairwise cutoff bond detection algorithm.

    Failure if number of detected clusters does not match known composition of given structure.

    Args:
        cluster: isolated structure to run TCC
    """
    try: x = xyz.read('clusters/%s.xyz' % cluster)
    except: return # temporary until we have sample coordinates for all known structures!

    found = run_static_tcc_voronoi(x)['Number']
    for component,n in structures.composition[cluster].items():
        assert n == int(found[component])
