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

	* 	CMake makefile generator
	*	Tetgen library 
	*	Armadillo
	*	VTK development library (for VTK import and export)
	* 	hdf5-dev (for importing Dream3D files)
	* 	libproj-dev (required by VTK)

To install all dependencies, simply run

	sudo apt-get install cmake libhdf5-dev libarmadillo-dev libvtk6-dev libtet1.5-dev libproj-dev

from the command line. 

To compile Voxel2Tet, first download Voxel2Tet source code. Using git, this is done by

$ git clone ....

Then create a build directory, e.g.

$ mkdir -p ~/bin/build/Voxel2Tet

Run CMake

$ cmake ~/dev/Voxel2Tet

For the standard "release" compilation, no extra flags are needed. The executable files are located in the Examples subdirectory.

Windows
-------

To do...

Mac OS
------

To do...