add_library(Voxel2TetLib
        Surface.cpp
        Voxel2Tet.cpp
        PhaseEdge.cpp
        Volume.cpp
        )

target_include_directories(Voxel2TetLib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/Export
        ${CMAKE_CURRENT_SOURCE_DIR}/General
        ${CMAKE_CURRENT_SOURCE_DIR}/Import
        ${CMAKE_CURRENT_SOURCE_DIR}/Mesh
        ${CMAKE_CURRENT_SOURCE_DIR}/Pipeline
        )

if (VOXEL2TET_EXTERNAL_INCLUDE_DIRS)
    target_include_directories(Voxel2TetLib SYSTEM PUBLIC ${VOXEL2TET_EXTERNAL_INCLUDE_DIRS})
endif ()

if (VOXEL2TET_EXTERNAL_LIBS)
    target_link_libraries(Voxel2TetLib PUBLIC ${VOXEL2TET_EXTERNAL_LIBS})
endif ()

target_compile_options(Voxel2TetLib PRIVATE -Wall -Wno-comment)

set(subdirs Export General Import Mesh Pipeline)
foreach (subdir ${subdirs})
    add_subdirectory(${subdir})
endforeach (subdir)
