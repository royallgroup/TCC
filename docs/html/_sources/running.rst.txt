Running the TCC
=================

Running the TCC requires 4 files, the binary, inputparameters.ini, box.txt and an input xyz file.

* A binary called ``tcc`` will be generated in the bin folder.
* inputparameters.ini sets some required parameters for the analysis, an example is given in the samples folder. Details of the parameters are given on the ::inputparameters page.
* box.txt gives the boundaries of the simulation box for systems with periodic boundaries. An example can be found in the examples folder. The specification of the box file is given below.
* The input configuration should be given as an XYZ file, it is not required that the XYZ file have the same number of particles in each frame. The TCC is limited to reading in a maximum of 1000 xyz frames.

Specifying periodic boundary conditions
----------------------------------------

* If the system does not have periodic boundary conditions set PBCs to 0 in inputparamers.in
* For any system with periodic boundaries the boundaries are stored in the box file specified in inputparamers.in. The first line of the box file should be a title line (usually column headings) which will be ignored. The specification of coordinates is detailed below for each system type.
* For NVT systems set the box_type parameter in inputparamers.in to 1. The box size will be read at TCC initialisation from the second line of the box.txt file (first line is a comment). The syntax is "timestep sidex sidey sidez" where timestep = 0 and sidex sidey and sidez are the x y and z box side lengths.
* For NPT systems set the box_type parameter in inputparamers.in to 2. The box size will be read each xyz frame from the box.txt file. The syntax is "timestep sidex sidey sidez" with each timestep on a new line. There must be at least as many timesteps as frames in the xyz file.
* For systems with triclinic boundary conditions with tilt, set the box_type parameter in inputparamers.in to 3. The box size will be read each xyz frame from the box.txt file. The syntax is "timestep sidex sidey sidez tilt" where the tilt has a sign.

Excluding clusters from analysis
***********************************

Sometimes it is not necessary to detect all possible cluster types, some clusters may be more interesting than others. For example in the hard sphere system 10B and 13A are of particular interest since they are locally favoured structures. Selecting only a subset of clusters can significantly increase the speed of the TCC.

To analyse only a subset of the clusters set the **analyse_all_clusters** parameter in the ``inputparamters.ini`` file to 0. In the directory the where the TCC will execute, place the ``clusters_to_analyse.ini`` file - an example can be found in the examples folder. It is only necessary to specify the clusters you want to detect, any prerequisite smaller clusters will be automatically selected.