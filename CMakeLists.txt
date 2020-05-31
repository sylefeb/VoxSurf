cmake_minimum_required(VERSION 3.10)
project(pyvoxsurf)

add_subdirectory(lib/LibSL-small)
set(XTL_INCLUDE_DIRS lib/xtl/include lib/xtensor/include lib/xtensor-blas/include lib/xtensor-python/include)

find_package(HDF5 COMPONENTS CXX)
if (${HDF5_FOUND})
  add_executable(pyvoxsurf main.cpp)
  target_compile_features(pyvoxsurf PRIVATE cxx_std_17)
  target_include_directories(pyvoxsurf PRIVATE ${HDF5_INCLUDE_DIRS} ${XTL_INCLUDE_DIRS})
  target_link_libraries(pyvoxsurf PRIVATE voxelizer ${HDF5_LIBRARIES})
endif(${HDF5_FOUND})

add_library(voxelizer STATIC src/pyvoxsurf/voxelizer.cpp)
target_compile_features(voxelizer PRIVATE cxx_std_17)
target_link_libraries(voxelizer PUBLIC LibSL)
target_include_directories(voxelizer PRIVATE ${XTL_INCLUDE_DIRS})

set_target_properties(LibSL PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(voxelizer PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

find_package(Python COMPONENTS Interpreter NumPy)
add_subdirectory(lib/pybind11)
pybind11_add_module(pyvoxsurf src/pyvoxsurf/pyvoxsurf.cpp)
target_compile_features(pyvoxsurf PUBLIC cxx_std_14)
target_include_directories(
  pyvoxsurf PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                  ${XTL_INCLUDE_DIRS} xtensor-python/include)
target_link_libraries(pyvoxsurf PUBLIC voxelizer)
install(TARGETS pyvoxsurf DESTINATION ${Python_SITELIB})
