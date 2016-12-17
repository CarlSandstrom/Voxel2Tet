Introduction
============

Voxel2Tet converts voxel representations to tetrahedral mesh with smooth interfaces. The Algorithm first generates a smooth surface using and then coarsens the fine mesh until no more coarsening is permitted. The coarse surface mesh is the used to create a tetrahedral mesh by calling the TetGen library. The purpose of the software is two-fold. First, smooth surfaces are desirable in many context, such as during investigations of interfaces between grains in micromechanics. Second, a reduction of degrees of freedom (i.e. vertices for describing the geometry) for decreasing the computational effort used in computational analysis. That being said, a decrease in degrees of freedom is not guaranteed and depends on the original voxel data and the geometry of the object. Typically a very coarse input yields more degrees of freedom. Then again, a coarse voxel description is ill-suited for computations anyway.

The code is in C++ and is written in a modular way such that new features, such as new import and export filters, are easy to implement. 

The original intent was to allow for voxel representations of microstructures to be performed, but the software is general and in no way bound to microstructures.

Installation
============

Currently, Voxel2Tet compiles using either `gcc` or `clang`. However, the coding is pretty standard and modifying the code to compile using other C++ compilers should not pose much of a problem. Note that some compiler dependent features are used, e.g. in logging routines, the compiler specific function `__FUNCTION__` is used. This is simply avoided by not compiling using logging (if you are not a developer, logging is not needed and it is turned off by default).

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

	$ git clone https://github.com/CarlSandstrom/Voxel2Tet.git

Then create a build directory, e.g.

	$ mkdir -p ~/bin/build/Voxel2Tet

Run CMake

	$ cmake ~/dev/Voxel2Tet/src/ && make

For the standard "release" compilation, no extra flags are needed. The executable files are located in the Examples subdirectory.

Fedora
------
The dependencies are the same as for Ubuntu. They can be installed by running

	$ sudo dnf install cmake hdf5-devel armadillo-devel vtk-devel tetgen-devel proj-devel

To compile, enter the cloned directory and run

	$ cmake ./src/ && make

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

Voxel2Tet can be used both as an API and as a command line driven tool. A few examples for API usage are available in the `./src/bin` folder.

Command line tool
-----------------

### Flags

After compiling Voxel2Tet, the command line tool is available in the `./bin/` directory. First, there are a few flags used to handle the input/output of the program as listed below.

Flag|Meaning
----|-------
-input _filename_    | _filename_ is the name of the input file. Currently supported files are `.dream3d` and `.vtu`.
-output _filename_  	| Here, _filename_ is the base filename for the output file without any extension. The extension is determined by the export flags (-export_XXXXX_).
-exportvtksurface  | Export the final surface in VTK format
-exportvtkvolume  	| Export the final volume in VTK format
-exportoff        	| Export the final surface in .OFF format
-exportoofem      	| Export final volume as input file for OOFEM
-exportabaqus     	| Export final volume as input file for Abaqus
-exportsteps      	| Export the result of each step in VTK format (mainly for debugging purposes)
-datacontainer _name_	| (Dream3D input) Name of data group, default 'VoxelDataContainer'. Note that this is used for compatibility with older versions of Dream3D.
-materialid _name_ 	| (Dream3D input) Field containing an identifier for the phase, default 'GrainIds'.  Note that this is used for compatibility with older versions of Dream3D.
-voxelcutout _arg_ | Only consider the subset of the input contained within the boundingbox defined by _arg_. Here, _arg_="[xmin ymin zmin xmax ymax zmax]" (include citations and brackets) where all data are integers.
-treatzeroasvoid   | Treats a material with ID 0 as void. By default, this is considered a solid.

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

The meaning of most of the following flags will only be discussed briefly and for full understanding, we refer to the paper by de Cougny.

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
-TOL_COL_MAXVOLUMECHANGE_FACTOR _value_ | Changes the default value of TOL_COL_MAXVOLUMECHANGE such that the value is _value_ times the size of a voxel. Default is 1.
-TOL_COL_MAXERROR_FACTOR _value_ | Changes the default value of TOL_COL_MAXERROR such that the value is _value_ times the size of a voxel. Default is 1.

A good choice of especially TOL_COL_MAXVOLUMECHANGE_FACTOR and TOL_COL_MAXERROR_FACTOR depends quite a lot of the size of the smallest objects. Choosing a large value gives a smoother surface and works fine if all objects consists of a much higher number of voxels. However, if the value is lower (or about the same) as the number of voxels of the smallest object, the can in effect vanish.

### Command line examples

Here follows a few examples of the use of the command line tool to get the user started.

The first simple example converts the input example `Subset.dream3d` to a volume with smoothed surfaces. From the `bin` subdirectory in the build directory, run the following line:

	$ ./Voxel2Tet -input Subset.dream3d -output /tmp/Default

Note that you have to add a path to the `Subset.dream3d` file. This produces three files in the `/tmp/` path: `Default.surface.vtp`, `Default.volume.vtu` and `Default.stat`. The names of the files pretty much says it all. The surfaces (inner and outer) are saved in the `Default.surface.vtp` file, the volumes are stored in `Default.volume.vtu` and statistics containing information on e.g. loss/gain in volume of the constituents are found in `Default.stat`. You can open the `.vtp` and `.vtu` files in Paraview for visualization. In paraview you can choose color by `Mat ID` the see the different materials.

Now, we will see what happens if we make the springs much stiffer. From the Examples subdirectory in the build directory, run the following line:

	$ ./Voxel2Tet -input Subset.dream3d -output /tmp/Stiff -spring_c_factor 0.1 -edge_spring_c_factor 0.1

Once again, you can open the resulting files (now called `Stiff.*` ) in Paraview. What you should see is a result not very different from voxel data (except that is it a tetrahedralization of voxelsurfaces).

Smoothing of very large structures can take several hours (or even days in some extreme cases), thus, it is of importance to have the correct flags set before running the program. Here, the flag voxelcutout comes in handy. The flag tells Voxel2Tet to only extract a small piece of the voxeldata and perform smoothing. This works well when working with microstructures since the geometry is somewhat equal throughout the volume. Run the following line:

	$ ./Voxel2Tet -input 40x40x40.dream3d -output /tmp/Cutout -voxelcutout "[10 10 10 20 25 30]"

So rather than smoothing a 40x40x40=64000 voxel file, we choose to smooth a 10*15*20=3000 voxel subset of the file and investigate the result. Now, we can alter the settings, run a test case and investigate the result without having to wait for the complete file to be converted.

API
---

In the `./src/Examples` subdirectory, there are a few examples on using the API. The main reason for using the API is to be able to create a tetrahedral mesh given an implicit geometry. By "Implicit" it is meant that an ID for a given coordinate is returned, hence, e.g. centre and radius for a sphere is not given directly to Voxel2Tet (which would be explicit).

To use the API, the function `LoadCallback` of a `Voxel2TetClass` type object is called with a function pointer and some basic information on the volume (origin, voxel size and the dimensions of the object) as arguments.

The API is probably best understood by carefully reading the source code in the Example subdirectory.

In the following, we will briefly discuss the code for `SingleSphere.cpp`. The file contains two functions: `GiveMaterialIDByCoordinateSphere` and the `main` function. `GiveMaterialIDByCoordinateSphere` is called from the Voxel2Tet library with a coordinate as arguments. Here, the function simply checks if the coordinate is located inside a sphere with its center in (0.5, 0.5, 0.5) and a radius of 0.25. If so, 1 is returned and otherwise a 2 is returned. This implies that the material inside the sphere will have ID 1 and the material outside ID 2. The IDs returned does not need to ordered in any way and is the same ID that is exported for it corresponding tetrahedral volume.

In the `main` function, we need to create an Options type object. This object takes care of command line arguments, so the compiled executable will be able to handle the same flags as discussed in the previous section. After the `Voxel2TetClass` object is created, the data is loaded by calling the `LoadCallback` member function. Then, the `Process` member function is called which performs the smoothing and coarsening of the surfaces. By calling the `Tetrahedralize` member function, the smoothed volume is tetrahedralized. What follows after that is simply exporting the result.

Some general remarks
--------------------

Some remarks on converting voxel data to tetrahedral meshes:

* Higher resolution of the voxel data implies better result.
* Stray voxels (single voxels of material A enclosed in other materials) gives a high density of tetrahedrons in that region.

Code documentation
=================
The code is documented using Doxygen. In order to generate the documentation, make sure that Doxygen is installed on you system. Then run doxygen using the configuration file `./src/doxygen.conf`. On linux, this is simply done using

	$ doxygen doxygen.conf

Coding conventions
==================

In order to maintain consistency in terms of coding convention, the tool "Uncrustify" (https://github.com/uncrustify/uncrustify) is used. If code is added or modified, run the tool using the file `uncrustify.cfg` in the project root directory as a configuration file. To run uncrustify on only one file, go to the root folder of Voxel2Tet and run the command (Linux):

	$ uncrustify -c ./uncrustify.cfg --replace [source file]

where [source file] is replaced by the source file (either .cpp or .h). In order to run uncrustify on the whole project, go to the project root folder and run

	$ find ./ -iname "*.cpp" -or -iname "*.h" | xargs -n 1 -I @ uncrustify -c ~/uncrustify.cfg --replace @

Thanks
======

Thanks to Tomas Akenine Möller for letting me use his code for triangle intersections (http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/).
