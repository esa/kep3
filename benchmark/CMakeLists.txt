# Find the xtensor bits.
find_package(Boost COMPONENTS program_options REQUIRED)

# And pagmo
find_package(Pagmo REQUIRED)

function(ADD_kep3_BENCHMARK arg1)
  add_executable(${arg1} ${arg1}.cpp)
  target_link_libraries(${arg1} PRIVATE kep3 Boost::boost Boost::program_options Pagmo::pagmo)
  target_compile_definitions(${arg1} PRIVATE XTENSOR_USE_FLENS_BLAS PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${kep3_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${kep3_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${kep3_CXX_FLAGS_RELEASE}>"
  )
  # Setup the C++ standard.
  target_compile_features(${arg1} PRIVATE cxx_std_20)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
endfunction()

ADD_kep3_BENCHMARK(convert_anomalies_benchmark)
ADD_kep3_BENCHMARK(propagate_lagrangian_benchmark)
ADD_kep3_BENCHMARK(lambert_problem_benchmark)
ADD_kep3_BENCHMARK(stm_benchmark)
ADD_kep3_BENCHMARK(leg_sims_flanagan_benchmark)
ADD_kep3_BENCHMARK(leg_sims_flanagan_hf_benchmark)
ADD_kep3_BENCHMARK(leg_sf_benchmark_simple)


