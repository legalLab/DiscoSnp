################################################################################
#  TOOLS GENERATION 
################################################################################
ADD_SUBDIRECTORY(kissnp2)
ADD_SUBDIRECTORY(kissreads2)
ADD_SUBDIRECTORY(QHC)
ADD_SUBDIRECTORY(read_file_names)
ADD_SUBDIRECTORY(create_coverage_h5_file)
# disco_haplotypes (C++ version of scripts/disco_haplotypes.py) needs a C++17 compiler (GCC >= 7) and
# CMake >= 3.8: otherwise it is skipped, and the pipelines use the python script (same outputs)
include (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG ("-std=c++17" COMPILER_SUPPORTS_CXX17)
if (COMPILER_SUPPORTS_CXX17 AND NOT CMAKE_VERSION VERSION_LESS 3.8)
    ADD_SUBDIRECTORY(disco_haplotypes)
else ()
    message (WARNING "disco_haplotypes not built (needs a C++17 compiler and CMake >= 3.8): run_discoSnp++.sh and run_discoSnpRad.sh -H will use scripts/disco_haplotypes.py (python >= 3.6 with numpy), which gives the same outputs, more slowly")
endif ()
