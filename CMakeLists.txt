cmake_minimum_required(VERSION 2.6)
project(VoxSurf)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/LibSL-small
  ${PROJECT_SOURCE_DIR}/LibSL-small/LibSL  
)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
  )

ADD_EXECUTABLE(VoxSurf
  main.cpp
  LibSL-small/LibSL/Math/Math.cpp
  LibSL-small/LibSL/Math/Vertex.cpp
  LibSL-small/LibSL/Mesh/Mesh.cpp
  LibSL-small/LibSL/Mesh/MeshFormat_stl.cpp
  LibSL-small/LibSL/Mesh/VertexFormat_dynamic.cpp
  LibSL-small/LibSL/System/System.cpp
  LibSL-small/LibSL/CppHelpers/CppHelpers.cpp
)
 
if(WIN32)
target_link_libraries(VoxSurf shlwapi)
endif(WIN32)
