# README #

Up-to-date version of the Topological Cluster Classification code.

### Compile ###
Enter the directory
    `src` 
and type 
    `gcc -O3 *.c -o TCC`

### Features ###

* Cubic boxes
* Orthogonal non-cubic boxes
* Triclinic boundary conditions

### Turn on the options for triclinic boundary conditions###

* Set the BOX option to 3 (triclinic with tilt)
* Provide an input file containing a one-line header and one single line for every configuration.
* Every line contains 5 entries: the iteration, Lx, Ly,Lz, and the tilt (with sign).


### Dynamic TCC ###

To run the dynamic TCC:

* First run the normal TCC, outputting the .bonds file.
* Turn on the dynamic flag in the TCC, output .dyn files.
* Run the post-processing analysis. This outputs .lives files which have the cumulative pdfs.