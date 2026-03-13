set(PipelineLibList
        ${CMAKE_CURRENT_SOURCE_DIR}/InputLoaderPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SurfaceBuilderPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SmoothingPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/TetPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ExportPipeline.cpp
        )

set(Voxel2TetLibList ${Voxel2TetLibList} ${PipelineLibList} PARENT_SCOPE)
