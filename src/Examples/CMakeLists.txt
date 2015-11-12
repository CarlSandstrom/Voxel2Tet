set(files TestOptions TestImporter TestArmadillo TestOctree TestVTK TestMeshData)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
