# README #

Simple python script to post process TCC output to find net TCC clusters.
Francesco Turci - February 2016

Requires python 2/3 and NumPy.

To change the list of clusters considered, edit the priority_list dictionary. The clusters listed first will be those of highest priority in the net calculation, those listed last will be lowest priority.

The required command line arguments is the stub of the xyz file name. The script will search the same directory for TCC raw output files for the clusters specified in the priority list.

Two files are output - inputname.xyz_net and inputname.xyz_gross.