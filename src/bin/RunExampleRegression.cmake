if (NOT DEFINED EXAMPLE_BINARY OR EXAMPLE_BINARY STREQUAL "")
    message(FATAL_ERROR "EXAMPLE_BINARY must be provided")
endif ()

if (NOT DEFINED OUTPUT_DIR OR OUTPUT_DIR STREQUAL "")
    message(FATAL_ERROR "OUTPUT_DIR must be provided")
endif ()

if (NOT DEFINED EXPECTED_FILES OR EXPECTED_FILES STREQUAL "")
    message(FATAL_ERROR "EXPECTED_FILES must be provided")
endif ()

file(REMOVE_RECURSE "${OUTPUT_DIR}")
file(MAKE_DIRECTORY "${OUTPUT_DIR}")

set(EXAMPLE_ARGS_LIST "")
if (DEFINED EXAMPLE_ARGS AND NOT EXAMPLE_ARGS STREQUAL "")
    string(REPLACE "|" ";" EXAMPLE_ARGS_LIST "${EXAMPLE_ARGS}")
endif ()

execute_process(
    COMMAND "${EXAMPLE_BINARY}" ${EXAMPLE_ARGS_LIST}
    WORKING_DIRECTORY "${OUTPUT_DIR}"
    RESULT_VARIABLE EXAMPLE_RESULT
    OUTPUT_VARIABLE EXAMPLE_STDOUT
    ERROR_VARIABLE EXAMPLE_STDERR
)

if (NOT EXAMPLE_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Regression command failed (${EXAMPLE_RESULT})\n"
        "Binary: ${EXAMPLE_BINARY}\n"
        "Args: ${EXAMPLE_ARGS_LIST}\n"
        "Stdout:\n${EXAMPLE_STDOUT}\n"
        "Stderr:\n${EXAMPLE_STDERR}\n"
    )
endif ()

string(REPLACE "|" ";" EXPECTED_FILES_LIST "${EXPECTED_FILES}")
set(MISSING_FILES "")
foreach (EXPECTED_FILE IN LISTS EXPECTED_FILES_LIST)
    if (NOT EXISTS "${OUTPUT_DIR}/${EXPECTED_FILE}")
        list(APPEND MISSING_FILES "${EXPECTED_FILE}")
    endif ()
endforeach ()

if (MISSING_FILES)
    message(FATAL_ERROR
        "Regression output validation failed.\n"
        "Binary: ${EXAMPLE_BINARY}\n"
        "Missing files: ${MISSING_FILES}\n"
        "Output directory: ${OUTPUT_DIR}\n"
    )
endif ()
