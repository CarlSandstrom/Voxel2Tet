set(files Cube MultiSphere SingleSphere Voxel2Tet)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
