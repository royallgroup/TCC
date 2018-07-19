Python Scripts
****************

A number of Python scripts are located in the ``tcc_python_scripts`` folder in the root directory of the TCC. These are provided for convinience interfacing with the TCC and post processing results. To in order to import these python scripts as modules it is necessary to add the ``tcc_python_scripts`` folder to your $PYTHONPATH environment variable. This can either be done by permenantly adding the path to the environment variables on your system or by dynamically adding the path each time the scripts are run. An example of the latter is::

    import sys
    import os
    sys.path.append(os.path.abspath("../../"))
    from tcc_scripts.file_readers import xyz
    
where you append the path of the main TCC directory.

Coordinate File Readers
==========================

The scripts in the file readers folder are designed to privde a unified interface for reading from different coordinate file types. These allow reading in configurations from XYZ, dynamo and .atom (LAMMPS) files.

:class:`~tcc_python_scripts.file_readers.snapshot.Snapshot` is the base class from which other file readers are defined. We define a *snapshot* as a single configuration of particles at one time. Regradless of the file format read, the data is stored in a common :class:`~tcc_python_scripts.file_readers.snapshot.Snapshot` object. Multiple snapshots may be present in an XYZ or atom file and these can be read in and each will be stored in a seperate :class:`~tcc_python_scripts.file_readers.snapshot.Snapshot` object.

Example XYZ file reader script
--------------------------------

The :meth:`tcc_python_scripts.file_readers.xyz.read()` method is a generator which returns sequential :class:`~tcc_python_scripts.file_readers.snapshot` objects from a file. By default the genertor will iterate over all snapshots in the file returning them sequentually. :: 

    >>> for frame in xyz.read("test/integration_tests/basic_voronoi/sample.xyz"):
    ...     print(frame.num_particles)
    59
    57
    58

To load only some of the frames you can use the parameter num_frames which will load that many frames from the begining of the file. ::
    
    >>> for frame in xyz.read("test/integration_tests/basic_voronoi/sample.xyz", num_frames=1):
    ...     print(frame.num_particles)
    59

For more complex operations such as reading every other frame you should load all frames and only operate on those frames needed. ::

    >>> num_frames_loaded = 0
    ... for frame in xyz.read("test/integration_tests/basic_voronoi/sample.xyz"):
    ...     if num_frames_loaded %2 == 0:
    ...         print(frame.num_particles)
    ...     num_frames_loaded += 1
    59
    58
    
To retrieve all snapshots at once use the python :class:`list` function to return a list continaing all of the snapshots. ::

    >>> from file_readers import xyz
    ... 
    ... particle_coordinates = list(xyz.read("../../test/integration_tests/basic_voronoi/sample.xyz", num_frames=2))
    ... print(particle_coordinates)
    [<snapshot n=59 t=None>, <snapshot n=57 t=None>]
    
    
Python Wrapper
===============

The TCC Python wrapper is designed to be a lightweight way of automating simple TCC analyses. It can analyse a single configuration, this is most easily loaded using the above file readers.

The :meth:`tcc_python_scripts.tcc.wrapper.TCCWrapper.run()` method will run the TCC with the provided coordinates and return the cluster count after completion. By default the script will run the TCC in a temporary
directory and deleting it after the run is complete. To save the results, for example raw or cluster xyz files you can provide a directory to the :meth:`tcc_python_scripts.tcc.wrapper.TCCWrapper.run()`
method where the results will be saved.

Input parameters are set to default values unless specified. They can be set by adding values to the releveant input_parameters attribute. For details on in input parameters see: :doc:`tcc_input_parameters`.

    
Example Wrapper Script
------------------------
::

    from tcc import wrapper
    from file_readers import xyz
    
    # Open a TCCWrapper object - this holds information about the simulation we want to run
    TCC_setup = wrapper.TCCWrapper()
    
    # Get the box size. This can be read from a file or input manually
    box = [10, 10, 10]
    
    # Get the coordinates. The file_readers scripts are a good way to read in coordinates from a file.
    particle_coordinates = list(xyz.read("../../test/integration_tests/basic_voronoi/sample.xyz", num_frames=1))[0].particle_coordinates
    
    TCC_setup.input_parameters['Run']['Frames'] = 1
    TCC_setup.input_parameters['Run']['Frames'] = 1
    results = TCC_setup.run(box, particle_coordinates, silent=False)
    
    print(results['Number of clusters'])
    print(results['Mean Pop Per Frame'])
    
Net TCC
========

Simple python script to post process TCC output to find net TCC clusters. Writteb by Francesco Turci - February 2016, edited by Peter Crowther - March 2018.

Description of net clusters
---------------------------------

The cluster populations output by the TCC are gross populations. This means that every particle in each cluster is reported. Some of these results are not very interesting since we know that the incidence of smaller clusters is much higher than that of larger clusters, e.g. almost everything is always in a 5A, 6A and 7A cluster.

A different measure of clusters is to consider the largest cluster each particle occurs in. If a particle is in an 11A and a 5A then we report only the 11A. We call these the net cluster populations.

This definition relies on a hierarchy of cluster types which determines which is the "most important" cluster a particle can appear in. Usually we define this priority list as the lowest energy structure for each number of particles in decreasing order of particle size for the system being considered. We give some examples of such priority lists for common structures below.

Using the net cluster script
----------------------------------
Requires Python 3 and NumPy (Python 2 may be supported but is untested).

The list of clusters considered is determined by the priority list. The clusters listed first will be those of highest priority in the net calculation, those listed last will be lowest priority.

The code requires a TCC raw file for each cluster specified in the priority list. The net script can be run directly from the command line or by
calling the :func:`~tcc_python_scripts.net_clusters.net.net_cluster_calculation` function.

To run from the command line
-------------------------------

The required command line argument is the directory containing the raw files and the priority list. For example::

    ./net.py ./raw_output  (FCC, 13A, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c)

The priority list must have the cluster names spelled exactly as the extensions on the raw files and the list must be in round brackets. The results are averaged over all frames in the input files and are output as a text file.

To run from a Python script
-----------------------------

::

    from from 
    from tcc_python_scripts.net_clusters import net
    net.net_cluster_calculation("./raw_output, [FCC, 13A, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c])
    

Cluster priority lists for common systems
-------------------------------------------

Hard Spheres: (FCC, 13A, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c)

Kob-Andersen: (13K, 12K, 11A, 10K, 9K, 8K, 7K, sp4c, sp3c)

Wahnstrom: (13A, 12B, 11W, 10B, 9B, 8A, sp5c, sp4c, sp3c)

Lennard-Jones: (13A, 12B, 11C, 10B, 9B, 8B, sp5c, sp4c, sp3c)

Sticky Spheres: (FCC, HCP, 13B, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c)