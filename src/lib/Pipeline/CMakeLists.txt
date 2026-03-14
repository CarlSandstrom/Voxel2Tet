target_sources(Voxel2TetLib
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/InputLoaderPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SurfaceBuilderPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SmoothingPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/TetPipeline.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ExportPipeline.cpp
        )
