cmake_minimum_required(VERSION 3.14)

project(scran_aggregate
    VERSION 0.2.3
    DESCRIPTION "Aggregate expression values across cells"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(scran_aggregate INTERFACE)
add_library(libscran::scran_aggregate ALIAS scran_aggregate)

target_include_directories(scran_aggregate INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scran_aggregate>)
target_compile_features(scran_aggregate INTERFACE cxx_std_17)

# Dependencies
option(SCRAN_AGGREGATE_FETCH_EXTERN "Automatically fetch scran_aggregate's external dependencies." ON)
if(SCRAN_AGGREGATE_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 3.0.0 CONFIG REQUIRED)
    find_package(tatami_tatami_stats 1.1.0 CONFIG REQUIRED)
endif()

target_link_libraries(scran_aggregate INTERFACE tatami::tatami tatami::tatami_stats)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(SCRAN_AGGREGATE_TESTS "Build scran_aggregate's test suite." ON)
else()
    option(SCRAN_AGGREGATE_TESTS "Build scran_aggregate's test suite." OFF)
endif()

if(SCRAN_AGGREGATE_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/scran_aggregate)

install(TARGETS scran_aggregate
    EXPORT scran_aggregateTargets)

install(EXPORT scran_aggregateTargets
    FILE libscran_scran_aggregateTargets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_aggregate)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_aggregateConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_aggregate)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_aggregateConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_aggregateConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_scran_aggregateConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_scran_aggregate)
