set(files TestOptions TestImporter TestOctree TestVTK TestMeshData TestVoxel2Tet TestCallbackImporter TestCube TestOpenMP)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
