# Configure the version file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/_version.py" @ONLY)   

# The list of pydsyre's Python files.
set(kep3_PYTHON_FILES __init__.py _version.py test.py)

# Copy the python files in the current binary dir,
# so that we can import pydesyre from the build dir.
# NOTE: importing from the build dir will work
# only on single-configuration generators.
foreach(kep3_PYTHON_FILE ${kep3_PYTHON_FILES})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${kep3_PYTHON_FILE}"
        "${CMAKE_CURRENT_BINARY_DIR}/${kep3_PYTHON_FILE}" COPYONLY)
endforeach()

# Core module.
Python3_add_library(core MODULE WITH_SOABI
    core.cpp
    docstrings.cpp
)

target_link_libraries(core PRIVATE kep3::kep3)
target_link_libraries(core PRIVATE "${pybind11_LIBRARIES}")

target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}" "${Python3_INCLUDE_DIRS}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")
target_compile_options(core PRIVATE
    "$<$<CONFIG:Debug>:${KEP3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${KEP3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${KEP3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${KEP3_CXX_FLAGS_RELEASE}>"
    
)
target_include_directories(core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
set_target_properties(core PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(core PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
target_compile_features(core PRIVATE cxx_std_17)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Installation setup.
if(PYKEP_INSTALL_PATH STREQUAL "")
    message(STATUS "pykep will be installed in the default location: ${Python3_SITEARCH}")
    set(_PYKEP_INSTALL_DIR "${Python3_SITEARCH}/pykep")
else()
    message(STATUS "pykep will be installed in the custom location: ${PYKEP_INSTALL_PATH}")
    set(_PYKEP_INSTALL_DIR "${PYKEP_INSTALL_PATH}/pykep")
endif()

# Install the core module.
install(TARGETS core
    RUNTIME DESTINATION ${_PYKEP_INSTALL_DIR}
    LIBRARY DESTINATION ${_PYKEP_INSTALL_DIR}
)

# Install the Python files.
install(FILES ${PYKEP_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR})

unset(_PYKEP_INSTALL_DIR)