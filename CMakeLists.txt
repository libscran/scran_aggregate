cmake_minimum_required(VERSION 3.14)

project(aggregate_across_cells
    VERSION 1.0.0
    DESCRIPTION "Aggregate statistics across cells"
    LANGUAGES CXX)

add_library(aggregate_across_cells INTERFACE)
add_library(scran::aggregate_across_cells ALIAS aggregate_across_cells)

include(GNUInstallDirs)
target_include_directories(aggregate_across_cells INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/scran_aggregate_across_cells>)

target_compile_features(aggregate_across_cells INTERFACE cxx_std_17)

# Dependencies
option(AGGREGATE_ACROSS_CELLS_FETCH_EXTERN "Automatically fetch aggregate_across_cells's external dependencies." ON)
if(AGGREGATE_ACROSS_CELLS_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 3.0.0 CONFIG REQUIRED)
endif()

target_link_libraries(aggregate_across_cells INTERFACE tatami::tatami)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(AGGREGATE_ACROSS_CELLS_TESTS "Build aggregate_across_cells's test suite." ON)
else()
    option(AGGREGATE_ACROSS_CELLS_TESTS "Build aggregate_across_cells's test suite." OFF)
endif()
if(AGGREGATE_ACROSS_CELLS_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/scran_aggregate_across_cells)

install(TARGETS aggregate_across_cells
    EXPORT aggregate_across_cellsTargets)

install(EXPORT aggregate_across_cellsTargets
    FILE scran_aggregate_across_cellsTargets.cmake
    NAMESPACE scran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scran_aggregate_across_cells)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/scran_aggregate_across_cellsConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scran_aggregate_across_cells)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/scran_aggregate_across_cellsConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/scran_aggregate_across_cellsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/scran_aggregate_across_cellsConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scran_aggregate_across_cells)
