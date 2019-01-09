Setup 
======
The latest source code for the TCC can be downloaded from GitHub: https://github.com/royallgroup/TCC

The TCC uses the free open source build system **CMake** (https://cmake.org/) to generate makefiles.

To build on systems other than Windows
----------------------------------------

From a terminal in the TCC directory execute::

    mkdir build
    cd build
    cmake ..
    make
    make install

To build on Windows with MinGW
---------------------------------

From a command prompt in the TCC directory execute::
    
    mkdir build
    cd build
    cmake .. -G "MinGW Makefiles"
    mingw32-make.exe
    mingw32-make.exe install
    
Compilation with Visual Studio is not tested.

Specifying install directory
------------------------------

By default the compiled TCC binary is installed to the bin directory in the loacation of the source. To specify an alternative install directory for the binary, use the command::

    cmake -DCMAKE_INSTALL_PREFIX:PATH=<PATH>
    
where ``<PATH>`` is the desired install directory.

To build with debug symbols
----------------------------

Building with debug symbols may be useful for debugging or profiling. Compile using::
    
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    
Testing the TCC
-----------------

It is recommended you test the TCC to check it is compiled correctly on your system, this requires Python and the pytest, NumPy and Pandas libraries. The test will run a short configurration to check clusters are correctly detected.

Before running the tests it required to first install the tcc_python_scripts module to your Python envrionment as described in :ref:`installation`

Once you have built the TCC, cd to the main TCC directory and type ``pytest``