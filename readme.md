Introduction
============

Voxel2Tet converts voxel representations to tetrahedral mesh with smooth interfaces. The Algorithm first generates a smooth surface using and then coarsens the fine mesh until no more coarsening is permitted. The coarse surface mesh is the used to create a tetrahedral mesh by calling the TetGen library. 

The code is in C++ and is written in a modular way such that new featers, such as new import and export filters, are easy to implement. It is text-based and should thus work on all common platforms.

The original intent was to allow for voxel representations of microstructures to be performed, but the software is general and in no way bound to microstructures.

Installation
============

Ubuntu
------

Requirements:

 - CMake makefile generator
 - Tetgen library 
 - Armadillo
 - VTK development library (for VTK import and export)
 - hdf5-dev (for importing Dream3D files)
 - libproj-dev (required by VTK)

To install all dependencies, simply run

	sudo apt-get install cmake libhdf5-dev libarmadillo-dev libvtk6-dev libtet1.5-dev libproj-dev

from the command line. 

To compile Voxel2Tet, first download Voxel2Tet source code. Using git, this is done by

	$ git clone git@github.com:CarlSandstrom/Voxel2Tet.git

Then create a build directory, e.g.

	$ mkdir -p ~/bin/build/Voxel2Tet

Run CMake

	$ cmake ~/dev/Voxel2Tet/src/

For the standard "release" compilation, no extra flags are needed. The executable files are located in the Examples subdirectory.

Windows
-------

To do...

Mac OS
------

To do...

Files
=====

./uncrustify.cfg 	Configuration file for uncrustify

./src/	Root location for all source files
./src/lib 	Runtime library
./src/bin	Source code for executables. 

./src/bin/Voxel2Tet.cpp 	Source code for main executable

./src/bin/Examples/SingleSphere.cpp
./src/bin/Examples/RandomSpheres.cpp

./ExampleInput/	Example files containing voxel data. Test main executable with these files.

Code documentation
=================
The code is documented using Doxygen. In order to generate the documentation, make sure that Doxygen is installed on you system. Then run doxygen using the configuration file ./src/doxygen.conf. On linux, this is simply done using

	$ doxygen doxygen.conf


Coding conventions
==================

In order to maintain consistency in terms of coding convention, the tool "Uncrustify" (https://github.com/uncrustify/uncrustify) is used. If code is added or modified, run the tool using the file "uncrustify.cfg" in the project root directory as a configuration file. To run uncrustify on only one file, go to the root folder of Voxel2Tet and run the command (Linux):

	$ uncrustify -c ./uncrustify.cfg --replace [source file]

where [source file] is replaced by the source file (either .cpp or .h). In order to run uncrustify on the whole project, go to the project root folder and run

	$ find ./ -iname "*.cpp" -or -iname "*.h" | xargs -n 1 -I @ uncrustify -c ~/uncrustify.cfg --replace @

Thanks
======

Thanks to Tomas Akenine Möller for letting me use his code (http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/).

