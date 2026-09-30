project(disco_haplotypes)

################################################################################
# TOOL: C++ implementation of scripts/disco_haplotypes.py (no GATB dependency)
################################################################################

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

file (GLOB_RECURSE  ProjectFiles  src/*.cpp)
add_executable (${PROJECT_NAME}  ${ProjectFiles})

# the same floating point results as the python script (no fused multiply-add)
target_compile_options (${PROJECT_NAME} PRIVATE -ffp-contract=off)

find_package (OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries (${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif ()

INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)
