set(files Cube MultiSphere SingleSphere FiberousMaterial Voxel2Tet)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
