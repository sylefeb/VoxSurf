cmake_minimum_required(VERSION 2.6)
project(VoxSurf)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/LibSL-small
  ${PROJECT_SOURCE_DIR}/LibSL-small/src/
  ${PROJECT_SOURCE_DIR}/LibSL-small/src/LibSL  
)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
  )

ADD_EXECUTABLE(VoxSurf
  main.cpp
  LibSL-small/src/LibSL/Math/Math.cpp
  LibSL-small/src/LibSL/Math/Vertex.cpp
  LibSL-small/src/LibSL/Mesh/Mesh.cpp
  LibSL-small/src/LibSL/Mesh/MeshFormat_stl.cpp
  LibSL-small/src/LibSL/Mesh/VertexFormat_dynamic.cpp
  LibSL-small/src/LibSL/System/System.cpp
  LibSL-small/src/LibSL/CppHelpers/CppHelpers.cpp
)
 
if(WIN32)
target_link_libraries(VoxSurf shlwapi)
endif(WIN32)
