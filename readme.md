Introduction
============

Voxel2Tet converts voxel representations to tetrahedral mesh with smooth interfaces. The Algorithm first generates a smooth surface using and then coarsens the fine mesh until no more coarsening is permitted. The coarse surface mesh is the used to create a tetrahedral mesh by calling the TetGen library. The purpose of the software is two-fold. First, smooth surfaces are desirable in many context, such as during investigations of interfaces between grains in micromechanics. Second, a reduction of degrees of freedom (i.e. vertices for describing the geometry) for decresing the computational effort used in computational analysis. That being said, a decrease in degrees of freedom is not guaranteed and depends on the original voxel data and the geometry of the object. Typically a very coarse input yields more degrees of freedom. Then again, a coarse voxel description is ill-suited for computations anyway.

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

File/path|Description
---------|-----------
./uncrustify.cfg |	Configuration file for uncrustify
./src/ |	Root location for all source files
./src/lib |	Runtime library
./src/bin	| Source code for executables.
./src/bin/Voxel2Tet.cpp |	Source code for main executable file.
./src/bin/SingleSphere.cpp |	Example of a single sphere contained in a cube.
./src/bin/MultiSphere.cpp |	Example of several, randomly places spheres in a block.
./src/bin/FiberousMaterial.cpp |	Example of a fiberous material.
./ExampleInput/ |	Example files containing voxel data. Test main executable with these files.

Usage
=====

Voxel2Tet can be used both as an API and as a command line driven tool. A few examples for API usage is avalible in the ./src/bin folder.

Command line tool
-----------------

After compiling Voxel2Tet, the command line tool is available in the ./Examples/ directory. First, there are a few flags used to handle the input/output of the program as listed below.

Flag|Meaning
----|-------
-output _filename_  	| Here, _filename_ is the base filename for the output file without any extension. The extension is determined by the export flags (-export_XXXXX_).
-exportvtksurface  | Export the final surface in VTK format
-exportvtkvolume  	| Export the final volume in VTK format
-exportoff        	| Export the funal surface in .OFF format
-exportoofem      	| Export final volume as input file for OOFEM
-exportabaqus     	| Export final volume as input file for Abaqus
-exportsteps      	| Export the result of each step in VTK format (mainly for debugging purposes)
-datacontainer _name_	| (Dream3D input) Name of data group, default 'VoxelDataContainer'. Note that this is used for compatibility with older versions of Dream3D.
-materialid _name_ 	| (Dream3D input) Field containing an identifier for the phase, default 'GrainIds'.  Note that this is used for compatibility with older versions of Dream3D.
-voxelcutout _arg_ | Only consider the subset of the input contained within the boundingbox defined by _arg_. Here, _arg_="[xmin ymin zmin xmax ymax zmax]" (include citations and brackets) where all data are integers.

Some more advanced flags for determining the behavior the smoothening algorithm are also available. For clarity, we first want to inform the reader that the smoothening algorithm consists of two parts. The first part is the smoothening part, where all vertices are moved in order to produce a smooth surface. Note that here all vertices are preserved. The second part is the mesh coarsening part where triangles are collapsed in order to reduce the number of vertices used and to smooth the surface further.

Before smoothing any surfaces, Voxel2Tet identifies edges. Here, an edge is where three or more materials meet. Flags for determining the behavior of the first part are shown in the table below.

Flag|Meaning
----|-------
-spring_c _value_           | Determines how far a vertex can be displaced before penalized. Default depends on the size of the voxels and is computed such that  a displacement is penalized when a vertex is moved farther away from its original position than the characteristic length of a voxel.
-spring_c_factor _value_    | Changes the value of spring_c to fit a multiple _value_ of the characteristic length of a voxel. Default is 1.
-spring_alpha _value_       | Determines penalization for vertices moving too far from the original position. Default is 2.
-edge_spring_c _value_      | See spring_c flag.
-edge_spring_factor _value_ | See spring_factor flag. Default is 1.
-edge_spring_alpha _value_  | See spring_alpha flag. Default is 3.

The mesh coarsening part is a variation of the mesh coarsening algorithm proposed by H.L. de Cougny (1998). Some features has been added to the de Cougny algorithm: 1) The new algorithm allow for edges 2) For each collapse of a vertex, the error (change of volume) is computed. A collapse implying a to large change in volume is not performed. 3) The error from a collapse is associated with the vertices affected by the collapse. This error is propagated for each collapse and there is a maximum threshold in the accumulated errors that cannot be exceeded.

The meaning of these flags will only be discussed briefly and for full understanding, we refer to the paper by de Cougny.

Flag|meaning
----|-------
-TOL_FLIP_MAXAREACHANGE _value_ | Largest change in area due to a flip of a shared edge. Default is 1e-2.
-TOL_FLIP_SMALLESTAREA _value_ | Smallest allowed area of a triangle after a flip. Default is 1e-8.
-TOL_FLIP_MAXNORMALCHANGE _value_ | Largest change in direction of normal of triangles due to a flip (radians). Default is 0.35.
-TOL_FLIP_MAXNORMALDIFFERENCE _value_ | Maximum angle of normals of two triangle to allow flipping (radians). Default is 0.175
-TOL_COL_SMALLESTAREA _value_ | Smallest allowed area of a triangle after collapsing. Default is 1e-8.
-TOL_COL_MINANGLE _value_ | Smallest inner angle of a triangle after collapsing (radians). Default is 0.35.
-TOL_COL_MAXNORMALCHANGE _value_ | Maximum change in normal direction due to a collapse (radians). Default is 0.35.
-TOL_COL_CHORD_MAXNORMALCHANGE _value_ | Maximum change in normal of a chord due to a collapse (radians). Default is 0.175.
-TOL_COL_MAXVOLUMECHANGE _value_ | Maximum change in volume due to one collapse. Default is 2 voxels.
-TOL_COL_MAXERROR _value_ | Largest accumulated error in a vertex. Default is 10 voxels.
-TOL_COL_MAXVOLUMECHANGE_FACTOR _value_ | Changes the default value of TOL_COL_MAXVOLUMECHANGE such that the value id _value_ times the size of a voxel. Default is 1.
-TOL_COL_MAXERROR_FACTOR _value_ | Changes the default value of TOL_COL_MAXERROR such that the value id _value_ times the size of a voxel. Default is 1.

### Command line examples

Here follows a few examples of the use of the command line tool to get the user started.

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

Thanks to Tomas Akenine Möller for letting me use his code for triangle intersections (http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/).
