cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(spline LANGUAGES CXX)

include(GNUInstallDirs)
set(PACKAGE_NAME spline)
set(USING_PYTHON3 True)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")

foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Library target
add_library(spline INTERFACE)
target_include_directories(spline
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(spline INTERFACE cxx_std_17)
target_link_libraries(spline INTERFACE)



install(TARGETS spline
  EXPORT spline-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/spline"
)

# Install targets file
install(EXPORT spline-targets
  FILE ${PACKAGE_NAME}Targets.cmake
  DESTINATION ${INSTALL_CMAKE_DIR}/cmake/spline
)

# Install splineConfig.cmake
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PACKAGE_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}/cmake/spline
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
  DESTINATION ${INSTALL_CMAKE_DIR}/cmake/spline
)

export(TARGETS spline
  FILE ${PACKAGE_NAME}Config.cmake)

# Install headers
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/spline.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)



set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE spline)
