cmake_minimum_required(VERSION 2.6)
project(voxsl)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/libsl-small
  ${PROJECT_SOURCE_DIR}/libsl-small/LibSL  
)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/path.h"
  )

ADD_EXECUTABLE(voxsl
  main.cpp
  libsl-small/LibSL/Math/Math.cpp
  libsl-small/LibSL/Math/Vertex.cpp
  libsl-small/LibSL/Mesh/Mesh.cpp
  libsl-small/LibSL/Mesh/MeshFormat_stl.cpp
  libsl-small/LibSL/Mesh/VertexFormat_dynamic.cpp
  libsl-small/LibSL/System/System.cpp
  libsl-small/LibSL/CppHelpers/CppHelpers.cpp
)
 
if(WIN32)
target_link_libraries(voxsl shlwapi)
endif(WIN32)
