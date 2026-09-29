set(files Cube MultiSphere SingleSphere FiberousMaterial Voxel2Tet)

foreach (file ${files})
    add_executable(${file} ${file}.cpp)
    target_link_libraries(${file} PRIVATE Voxel2TetLib)
endforeach (file)

if (BUILD_TESTING)
    set(REGRESSION_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/RunExampleRegression.cmake")
    set(REGRESSION_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/ctest_outputs")
    set(CUBE_EXPECTED_FILES "Cube.vtp|Cube.vtu")

    add_test(
        NAME regression_cube_example
        COMMAND ${CMAKE_COMMAND}
            -DEXAMPLE_BINARY=$<TARGET_FILE:Cube>
            -DOUTPUT_DIR=${REGRESSION_OUTPUT_DIR}/Cube
            "-DEXPECTED_FILES=${CUBE_EXPECTED_FILES}"
            -P ${REGRESSION_SCRIPT}
    )

    set(EXAMPLE_INPUT_16X16X16 "${CMAKE_SOURCE_DIR}/../ExampleInput/16x16x16.dream3d")
    if (EXISTS "${EXAMPLE_INPUT_16X16X16}")
        set(CLI_EXAMPLE_ARGS "-input|${EXAMPLE_INPUT_16X16X16}|-output|dream3d_regression")
        set(CLI_EXPECTED_FILES
            "dream3d_regression.surface.vtp|dream3d_regression.volume_0.vtu|dream3d_regression.stat")
        add_test(
            NAME regression_voxel2tet_cli_16x16x16
            COMMAND ${CMAKE_COMMAND}
                -DEXAMPLE_BINARY=$<TARGET_FILE:Voxel2Tet>
                -DOUTPUT_DIR=${REGRESSION_OUTPUT_DIR}/Voxel2Tet16x16x16
                "-DEXAMPLE_ARGS=${CLI_EXAMPLE_ARGS}"
                "-DEXPECTED_FILES=${CLI_EXPECTED_FILES}"
                -P ${REGRESSION_SCRIPT}
        )
        set_tests_properties(regression_voxel2tet_cli_16x16x16 PROPERTIES TIMEOUT 120)
    else ()
        message(WARNING "Skipping regression_voxel2tet_cli_16x16x16: missing ${EXAMPLE_INPUT_16X16X16}")
    endif ()
endif ()
