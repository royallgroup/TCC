TCC input parameters
===========================

Paramters for the TCC are set in the file inputparameters.ini. Examples can be found in the examples folder and the format in the example
files should be followed.

The input parameters are broken down into ``Box``, ``Run``, ``Simulation`` and ``Ouptut`` sections. For each section the parameters are given
along with the default values in **bold** and a brief description.

None of the parameters are mandatory - if a parameter is not provided then its default value will be used. 

Box Parameters
---------------

Specifies how to read simulation box

:box_type:			**1**				; 1 if NVT, 2 if NPT, 3 if triclinc with tilt

:box_name:			**box.txt**     	; name of parameters file for box size

Run Parameters
------------------
These settings describe the configuration or configurations to be loaded from the xyz input file.

:xyzfilename:         **sample.xyz**    ; File name of the xyz file to be analysed.

:frames:              **1**             ; Number of frames to load and analyse


Simulation
------------
These parameters control the interactions between the particles and depend on the type of system you are analysing. Most of these paramters
control the construction of the bond network which determines which particles are neighbours.

:rcutAA:				**1.8**	    ; maximum bond lengths between two A species particles, the cutoff is always applied whether Voronoi bonds are used or not.

:rcutAB:				**1.8**	    ; maximum bond lengths between an A and a B species particle, the cutoff is always applied whether Voronoi bonds are used or not.

:rcutBB:				**1.8**	    ; maximum bond lengths between two B species particles, the cutoff is always applied whether Voronoi bonds are used or not.

:min_cutAA:             **0.0**     ; minimum A-A bond length. Good for excluding overlapping particles in ideal gases.

:bond_type:			    **1**		; How bonds are detected. Set to 0 to use a simple bond length cutoff,  set to 1 to user Voronoi bond detection.

:PBCs:				    **1**       ; Set to 1 to use periodic boundary conditions when determining the bond network. Set to 0 to turn off perioid boundary conditions,

:voronoi_parameter:     **0.82**    ; Modified Voronoi Fc parameter. Parameter for controlling Voronoi bond detection.

:num_bonds:			    **50**	    ; Maximum number of bonds to one particle.

:cell_list:			    **0**		; Use cell list to calculate bond network - may speed up construction of bond network.

:analyse_all_clusters:  **1**       ; When set to 1 the TCC will analyse all possible clusters. If set to zero will read only the subset of clusters specified in ``clusters_to_analyse.ini``

Output	
---------
These parameters determine what files the TCC will output at the end of the run. The average populations of clusters for the whole trajectory will always be output.

:bonds: 				**0**		; Write bonds file which describes the bond network.

:clusts: 				**0**		; Write clusts files containing which list the ids of particles in all clusters

:raw: 				    **0**      	; Write raw cluster files. These specify whether a particle is or isn't in each structure.

:do_XYZ:                **0**       ; Write clusters to xyz files. This generates an XYZ file containing only those particles in each cluster type.

:11a: 				    **0**		; Write out an XYZ file containing only the central particles of 11A clusters.

:13a: 				    **0**		; Write out an XYZ file containing only the central particles of 13A clusters.

:pop_per_frame: 		**0**		; Write out the population of each cluster in each frame.