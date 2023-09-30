
add_library(kep3_test STATIC catch_main.cpp)
target_compile_options(kep3_test PRIVATE
  "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
)

# Setup the C++ standard.
target_compile_features(kep3_test PRIVATE cxx_std_20)
set_property(TARGET kep3_test PROPERTY CXX_EXTENSIONS NO)

function(ADD_kep3_TESTCASE arg1)
  add_executable(${arg1} ${arg1}.cpp)
  target_link_libraries(${arg1} PRIVATE kep3_test kep3 xtensor xtensor-blas Boost::boost)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
  )

  # Setup the C++ standard.
  target_compile_features(${arg1} PRIVATE cxx_std_20)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
  add_test(${arg1} ${arg1})
endfunction()

ADD_kep3_TESTCASE(convert_anomalies_test)
ADD_kep3_TESTCASE(epoch_test)
ADD_kep3_TESTCASE(planet_test)
ADD_kep3_TESTCASE(planet_keplerian_test)
ADD_kep3_TESTCASE(planet_jpl_lp_test)
ADD_kep3_TESTCASE(ic2par2ic_test)
ADD_kep3_TESTCASE(ic2eq2ic_test)
ADD_kep3_TESTCASE(eq2par2eq_test)
ADD_kep3_TESTCASE(propagate_lagrangian_test)
ADD_kep3_TESTCASE(propagate_keplerian_test)
ADD_kep3_TESTCASE(lambert_problem_test)