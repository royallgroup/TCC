# README #

Latest version of the Topological Cluster Classification code.

### Compilation ###
The TCC uses the free open source build system _cmake_ to generate makefiles.

##### To build on Windows with MinGW:

From a command prompt in the TCC directory execute
* `mkdir build`
* `cd build`
* `cmake .. -G "MinGW Makefiles`
* `mingw32-make.exe`

##### To build on other systems

From a terminal in the TCC directory execute
* `mkdir build`
* `cd build`
* `cmake ..`
* `make`

##### Testing the TCC

It is recommended you test the TCC to check it is compiled correctly on your system, this currently requires Python and the pytest module. The test will run a short configurration to check clusters are correctly detected.

Once you have built the TCC, cd to the test directory and type
'py.test test.py'

##### Running the TCC

Binaries will be generate in the bin folder.
Examples of the required input files are given in the examples folders.
To run the TCC the executable should be in the same directory as the input files.

### Setting up the simulation box ###

* For cubic boxes set the ISNOTCUBIC parameter in inputparamers.in to 0. The box size will be determined from the density given in the inputparamers.in file.
* For non-cubic NVT set the ISNOTCUBIC parameter in inputparamers.in to 1. The box size will be read at TCC initialisation from the second line of the box.txt file (first line is a comment). The syntax is "timestep sidex sidey sidez" where timestep = 0 and sidex y and z are the x y and z box side lengths.
* For non-cubic NPT set the ISNOTCUBIC parameter in inputparamers.in to 2. The box size will be read each timestep (each xyz frame) from the box.txt file. The syntax is "timestep sidex sidey sidez" with each timestep on a new line. There must be at least as many timesteps as frames in the xyz file.
* For Triclinic boundary conditions with tilt, set the ISNOTCUBIC parameter in inputparamers.in to 3. The box size will be read each timestep (each xyz frame) from the box.txt file. The syntax is "timestep sidex sidey sidez tilt" where the tilt has a sign.

### Dynamic TCC ###

The general schema to run the dynamic TCC is:

* First run the normal TCC, outputting the .bonds file.
* Turn on the dynamic flag in the TCC, output .dyn files.
* Run the post-processing analysis. This outputs .lives files which have the cumulative pdfs.
