# README #

Simple python script to post process TCC output to find net TCC clusters.
Francesco Turci - February 2016
Peter Crowther - March 2018

## Description of net clusters

The cluster populations output by the TCC are gross populations. This means that every particle in each cluster is reported. Some of these results are not very interesting since we know that the incidence of smaller clusters is much higher than that of larger clusters, e.g. almost everything is always in a 5A, 6A and 7A cluster.

A different measure of clusters is to consider the largest cluster each particle occurs in. If a particle is in an 11A and a 5A then we report only the 11A.

This definition relies on a hierarchy of cluster types which determines which is the "most important" cluster a particle can appear in. Usually we define this priority list as the lowest energy structure of each number of particles in decreasing order of particle size for the system being considered. We give some examples of such priority lists for common structures below.

## Using the net cluster script
Requires Python 3 and NumPy (Python 2 may be supported but is untested).

To change the list of clusters considered, edit the priority list dictionary. The clusters listed first will be those of highest priority in the net calculation, those listed last will be lowest priority.

The script must be run in the same directory as a folder called "raw_output" which contains the raw output files for each cluster type. The folder must contain a TCC raw file for each cluster specified in the priority list.

The required command line argument is the stub of the xyz file name.

The results are averaged over all frames in the input files and are output as a text file.

## Cluster priority lists for common systems

Hard Spheres: 'FCC', '13A', '12E', '11F', '10B', '9B', '8B', 'sp5c', 'sp4c', 'sp3c'
Kob-Andersen: '13K', '12K', '11A', '10K', '9K', '8K', '7K', 'sp4c', 'sp3c'
Wahnstrom: '13A', '12B', '11W', '10B', '9B', '8A', 'sp5c', 'sp4c', 'sp3c'
Lennard-Jones: '13A', '12B', '11C', '10B', '9B', '8B', 'sp5c', 'sp4c', 'sp3c'
Sticky Spheres: 'FCC', 'HCP', '13B', '12E', '11F', '10B', '9B', '8B', 'sp5c', 'sp4c', 'sp3c'

## Tests

The script comes with a test suite which requires the pytest library. This is located in /TCC/test/net_tcc_script.