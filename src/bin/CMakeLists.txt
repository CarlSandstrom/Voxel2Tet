set(files TestOctree TestMeshData Voxel2Tet TestCallbackImporter TestCube TestTriangleIntersectionOneShared)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
