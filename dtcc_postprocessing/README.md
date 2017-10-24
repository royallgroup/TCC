### README - DTCC Post Processing

The raw output from the dynamic TCC is not the final result. This post processing code needs to be run after the dynamic TCC in order to produce cumulative cluster lifetime distributions. This code serves to accumulate the lifetimes and also correct some biases present in the original code.

As for the TCC - iniparser must be built first before linking it to the code.

### Compile ###
Enter the directory
    `src/iniparser` 
and type 
	`make`
return to the src directory and type
    `gcc -O3 -lm *.c -o TCCpostprocessing -Liniparser/ -liniparser`

To run the dynamic TCC:

* Run the TCC once in static mode, outputting the .bonds file
* Run the TCC again in dynamic mode with the desired clusters set in dynamic memsize.dat. This should output .dyn files.
* Run this code, this will generate .lives files. These are the cumulative lifetime distributions.

#### Input files
This code requires the analysisparameters.ini file. If non-cubic boxes are used, it will also read in the box file specified in analysisparameters.ini. Analysis of systems with Triclinic boundary conditions is not currently supported.

The parameters in the analysisparameters.ini file are as follows:

* Analyse dynamic clusters - Set this to 1
* Maximum number of frames to remove a disappearance of a cluster - This allows the code to gloss over the temporary disappearance of clusters due to poor detection by the TCC or thermal fluctuations. If the code finds two clusters in the same place less than this number of frames apart then it assumes that they are one cluster and it has existed all along. Set this to ~ 1 tau alpha.
* remove instances where a subcluster becomes unbonded - This is described in Alex's thesis (pg. 200), in short, in the case above where a cluster disappears but then reappears within 1 tau alpha, this is a check that the cluster mostly stays together during the periods it is not detected by the TCC. Leave this on.
* re-write dyn_... file as a test - useful for debugging
* write processed clusters to new_... file - useful for debugging
* number of frames for N_clus^dyn/N to reach steady state - This accounts for the fact that clusters are not properly detected in the first X and last X frames of the trajectory due to the fact that you don't know what the clusters were doing before the start or end of the trajectory. Set this to the length of the longest lived cluster. If you dont have that much data you can decrease this but the smaller it is, the more it will skew the data towards sort lifetimes at the start and end.
* write cluster lifetime correlation functions - output the cumulative pdfs as _lives files.
* We don't know what the rest of the parameters do. Fiddling with them doesn't seem to make a great deal of difference so they probably aren't that important...