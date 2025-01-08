set(PYKEP_TRAJOPT_PYTHON_FILES __init__.py _direct_point2point.py _direct_pl2pl.py _mga.py _launchers.py)
install(FILES ${PYKEP_TRAJOPT_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/trajopt)

ADD_SUBDIRECTORY(gym)

