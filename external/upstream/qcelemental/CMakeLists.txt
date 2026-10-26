if(NOT (${CMAKE_DISABLE_FIND_PACKAGE_qcelemental}))
    include(FindPythonModule)
    find_python_module(qcelemental ATLEAST 0.50.4 QUIET)  # edit in codedeps
endif()

if(${qcelemental_FOUND})
    message(STATUS "${Cyan}Found qcelemental${ColourReset}: ${PY_qcelemental} (found version ${qcelemental_VERSION})")
    add_library(qcelemental_external INTERFACE)  # dummy

    if(${CMAKE_INSIST_FIND_PACKAGE_qcelemental})
        message(VERBOSE "Suitable qcelemental located externally as user insists")
    endif()

else()
    if(${CMAKE_INSIST_FIND_PACKAGE_qcelemental})
        message(FATAL_ERROR "Suitable qcelemental could not be externally located as user insists")
    endif()

    include(ExternalProject)
    message(STATUS "Suitable qcelemental could not be located, ${Magenta}Building qcelemental${ColourReset} instead.")

    file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _install_lib)

    set(_qcel_url
        https://github.com/MolSSI/QCElemental/archive/v0.50.4.tar.gz  # edit in codedeps
        )

    # qcelemental is pure Python, so installing it needs no pip, only its package dir plus the .dist-info
    #   that importlib.metadata reads for qcelemental.__version__. The source has no version in it
    #   (setuptools-scm takes it from git), so take it from the tag in the URL. For an untagged commit,
    #   set _qcel_version by hand instead (>= the find_python_module minimum above).
    string(REGEX REPLACE "^.*/v?([^/]+)\\.tar\\.gz$" "\\1" _qcel_version "${_qcel_url}")
    set(_qcel_distinfo ${CMAKE_CURRENT_BINARY_DIR}/qcelemental-${_qcel_version}.dist-info)
    file(WRITE ${_qcel_distinfo}/METADATA "Metadata-Version: 2.1\nName: qcelemental\nVersion: ${_qcel_version}\n")

    ExternalProject_Add(qcelemental_external
        URL ${_qcel_url}
        BUILD_ALWAYS 1
        CONFIGURE_COMMAND ""
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/qcelemental ${_install_lib}/qcelemental
                COMMAND ${CMAKE_COMMAND} -E copy_directory ${_qcel_distinfo} ${_install_lib}/qcelemental-${_qcel_version}.dist-info
        )
endif()
