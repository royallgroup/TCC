#XYZ file specification
The XYZ file format is a format to describe configurations of particle systems in 3D cartesian space.

This document sets out a specification for XYZ files to be read by the TCC, whether this be natively or through the python interface.

The XYZ file format does not have a strict specification which can make interpreting XYZ file programtically challenging. By specifying an exact format we hope to assist the user by defining the expected behavor of the program.

An XYZ file is comprised of "frames" which are a description of a particle system at a single point in time. Multiple frames may be appended in a file to describe a configuration over time.

## Structure of a frame
The first line of a frame specifies the number of particles (N) in the frame. It is an integer number. No other text is allowed on this line. 

The second line is a comment line. A comment may be placed here or the line may be left blank. This line is igored by the program.

There are then N lines, each of which describes the coordinates of a single particle. These lines consist of the identity of a particle followed by 3 spatial coordinates. No other text may be included in this line. 

The identity of a particle is specified by a single letter or number. The coordinates are given as floating point numbers. Each of these elements is separated by either a single space or single tab-space.

### Example of an XYZ frame
```
3
This is a comment line
A 5.67 -3.45 2.61
B 3.91 -1.91 4
A 3.2 1.2 -12.3
```

## Using multiple timesteps
If there are multiple timesteps then each timestep is appended directly after the last. It is not required that any quantities are conserved between timesteps (number of particles, particle identities etc.), each timestep is treated separately. It is not required to label or otherwise number frames although this is a good use of the comment line.

###Example of an XYZ file consisting of multiple timesteps
```
3
Frame 1
A 5.67 -3.45 2.61
B 3.91 -1.91 4
A 3.2 1.2 -12.3
4
Frame 2
B 5.47 -3.45 2.61
B 3.91 -1.93 3.1
A 3.2 1.2 -22.4
A 3.2 1.2 -12.3
3
Frame 3
1 5.67 -3.45 2.61
1 3.91 -1.91 4
2 3.2 1.2 -12.3
```