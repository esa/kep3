set(PYKEP_TRAJOPT_PYTHON_FILES 
    __init__.py 
    _direct_point2point.py 
    _direct_pl2pl.py 
    _mga.py _mga_1dsm.py
    _min_Bu_bu.py 
    _launchers.py 
    _mim.py
    _pl2pl_N_impulses.py
    _pontryagin_cartesian.py
    _primer_vector.py)
    
install(FILES ${PYKEP_TRAJOPT_PYTHON_FILES} DESTINATION ${_PYKEP_INSTALL_DIR}/trajopt)

ADD_SUBDIRECTORY(gym)