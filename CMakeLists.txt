cmake_minimum_required(VERSION 3.10)
project(VoxSurf)

add_subdirectory(LibSL-small)
set(XTL_INCLUDE_DIRS xtl/include xtensor/include xtensor-blas/include)

find_package(HDF5 COMPONENTS CXX)
if (${HDF5_FOUND})
  add_executable(VoxSurf main.cpp)
  target_compile_features(VoxSurf PRIVATE cxx_std_17)
  target_include_directories(VoxSurf PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(VoxSurf PRIVATE voxelizer ${HDF5_LIBRARIES})
endif(${HDF5_FOUND})

add_library(voxelizer STATIC voxelizer.cpp)
target_compile_features(voxelizer PRIVATE cxx_std_17)
target_link_libraries(voxelizer PUBLIC LibSL)
target_include_directories(voxelizer PRIVATE ${XTL_INCLUDE_DIRS})

set_target_properties(LibSL PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(voxelizer PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

find_package(Python COMPONENTS Interpreter NumPy)
add_subdirectory(pybind11)
pybind11_add_module(pyvoxel pyvoxel.cpp)
target_compile_features(pyvoxel PUBLIC cxx_std_14)
target_include_directories(
  pyvoxel PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS}
                  ${XTL_INCLUDE_DIRS} xtensor-python/include)
target_link_libraries(pyvoxel PUBLIC voxelizer)
install(TARGETS pyvoxel DESTINATION ${Python_SITELIB})
