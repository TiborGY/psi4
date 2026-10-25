# NumPy is needed at buildtime only for (1) the FindPython checks, (2) pybind11 version
#   selection, and (3) gau2grid code generation. Psi4 C++ doesn't use the NumPy C API.
#   So if NumPy is missing, build it from source and pip install it into the stage dir
#   alongside the psi4 pymodule, where buildtime, ctest, and runtime can all import it.
#
# Communicated to other projects:
# * psi4_NumPy_VERSION     -- version for pybind11 selection
# * psi4_NumPy_INCLUDE_DIR -- passed to psi4-core as Python_NumPy_INCLUDE_DIR hint
# * psi4_NumPy_STAGED      -- ON if NumPy lives in stage dir, so needs PYTHONPATH for buildtime import

file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

if(${Python_NumPy_FOUND} AND NOT (${CMAKE_DISABLE_FIND_PACKAGE_NumPy}))
    message(STATUS "${Cyan}Found NumPy${ColourReset}: ${Python_NumPy_INCLUDE_DIRS} (found version ${Python_NumPy_VERSION})")
    add_library(numpy_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_NumPy})
        message(VERBOSE "Suitable NumPy located externally as user insists")
    endif()

    set(psi4_NumPy_VERSION ${Python_NumPy_VERSION} CACHE INTERNAL "")
    # pass through any user hint; psi4-core otherwise detects NumPy itself
    set(psi4_NumPy_INCLUDE_DIR "${Python_NumPy_INCLUDE_DIR}" CACHE INTERNAL "")
    set(psi4_NumPy_STAGED OFF CACHE INTERNAL "")

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_NumPy})
        message(FATAL_ERROR "Suitable NumPy could not be externally located as user insists")
    endif()

    # catch NumPy built into stage by a previous configure of this build dir
    if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_NumPy}))
        include(FindPythonModule)
        find_python_module(numpy ATLEAST 2.0 QUIET)  # edit in codedeps
    endif()

    if(${numpy_FOUND})
        execute_process(COMMAND "${Python_EXECUTABLE}" "-c"
                                "import sys; \
                                 sys.path.insert(0, '${_install_lib}'); \
                                 import numpy; \
                                 print(numpy.get_include())"
            OUTPUT_VARIABLE _numpy_include
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "${Cyan}Found NumPy${ColourReset}: ${PY_numpy} (found version ${numpy_VERSION})")
        add_library(numpy_external INTERFACE)  # dummy

        set(psi4_NumPy_VERSION ${numpy_VERSION} CACHE INTERNAL "")
        set(psi4_NumPy_INCLUDE_DIR "${_numpy_include}" CACHE INTERNAL "")
        set(psi4_NumPy_STAGED ON CACHE INTERNAL "")

    else()
        if(WIN32)
            message(FATAL_ERROR "Suitable NumPy could not be located, and building NumPy is not supported on Windows. Install NumPy (e.g., `conda install numpy`)")
        endif()

        include(ExternalProject)
        message(STATUS "Suitable NumPy could not be located, ${Magenta}Building NumPy${ColourReset} instead.")

        # * Source build needs C & C++ compilers and Python headers (all already required by Psi4)
        #   and, through pip build isolation, PyPI access for meson-python & Cython.
        # * BLAS/LAPACK is NumPy's own detection (pkg-config, etc.), falling back to its slow lapack_lite.
        ExternalProject_Add(numpy_external
            DOWNLOAD_COMMAND ""
            CONFIGURE_COMMAND ""
            UPDATE_COMMAND ""
            BUILD_COMMAND ""
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
                            ${Python_EXECUTABLE} -m pip install numpy~=2.0  # edit in codedeps
                                                          --no-binary=numpy
                                                          --target=${_install_lib}
                                                          --no-deps
            USES_TERMINAL_INSTALL 1
            )

        get_filename_component(_numpy_include "${_install_lib}/numpy/_core/include" ABSOLUTE)
        set(psi4_NumPy_VERSION 2.0 CACHE INTERNAL "")
        set(psi4_NumPy_INCLUDE_DIR "${_numpy_include}" CACHE INTERNAL "")
        set(psi4_NumPy_STAGED ON CACHE INTERNAL "")
    endif()
endif()
