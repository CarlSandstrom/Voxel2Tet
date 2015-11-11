set(files TestOptions TestImporter TestArmadillo)

foreach(file ${files})
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} ${linkerlist})
endforeach(file)
