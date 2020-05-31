CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(LibSL-small)

SET(LibSL_DIR ${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(
			${CMAKE_CURRENT_SOURCE_DIR}/src
			${CMAKE_CURRENT_SOURCE_DIR}/src/LibSL
			${CMAKE_CURRENT_SOURCE_DIR}/src/LibSL/loki
			)

SET(LIBSL_SMALL_CORE
	src/LibSL/LibSL.h
	src/LibSL/CppHelpers/BasicParser.h
	src/LibSL/CppHelpers/CppHelpers.cpp
	src/LibSL/CppHelpers/CppHelpers.h
	src/LibSL/Errors/Errors.h
	src/LibSL/Geometry/AAB.h
	src/LibSL/Image/Filter.h
	src/LibSL/Image/Image.cpp
	src/LibSL/Image/Image.h
	src/LibSL/Image/ImageFormat_TGA.cpp
	src/LibSL/Image/ImageFormat_TGA.h
	src/LibSL/Image/tga.cpp
	src/LibSL/Image/tga.h
	src/LibSL/Math/Math.cpp
	src/LibSL/Math/Math.h
	src/LibSL/Math/Matrix4x4.h
	src/LibSL/Math/Tuple.h
	src/LibSL/Math/Vertex.cpp
	src/LibSL/Math/Vertex.h
	src/LibSL/Memory/Array.h
	src/LibSL/Memory/Array2D.h
	src/LibSL/Memory/Array3D.h
	src/LibSL/Memory/ArrayTools.h
	src/LibSL/Memory/Pointer.h
	src/LibSL/Mesh/Mesh.cpp
	src/LibSL/Mesh/Mesh.h
	src/LibSL/Mesh/MeshFormat_stl.cpp
	src/LibSL/Mesh/MeshFormat_stl.h
	src/LibSL/Mesh/VertexFormat.h
	src/LibSL/Mesh/VertexFormat_base.h
	src/LibSL/Mesh/VertexFormat_descriptor.h
	src/LibSL/Mesh/VertexFormat_dynamic.cpp
	src/LibSL/Mesh/VertexFormat_dynamic.h
	src/LibSL/Mesh/VertexFormat_reference.h
	src/LibSL/StlHelpers/StlHelpers.cpp
	src/LibSL/StlHelpers/StlHelpers.h
	src/LibSL/SvgHelpers/SvgHelpers.cpp
  src/LibSL/System/System.cpp
	src/LibSL/System/System.h
	src/LibSL/System/Types.h
	)

IF(NOT WIN32)
SET(LIBSL_SMALL_CORE
  ${LIBSL_SMALL_CORE}
	src/LibSL/Image/ImageFormat_JPG.cpp
	src/LibSL/Image/ImageFormat_JPG.h
	src/LibSL/Image/ImageFormat_PNG.cpp
	src/LibSL/Image/ImageFormat_PNG.h
)
ENDIF(NOT WIN32)

ADD_LIBRARY(LibSL STATIC ${LIBSL_SMALL_CORE})

IF(WIN32)
TARGET_LINK_LIBRARIES(LibSL shlwapi psapi)
ELSE(WIN32)
FIND_PACKAGE( ZLIB REQUIRED ) 
TARGET_LINK_LIBRARIES(LibSL png jpeg pthread ${ZLIB_LIBRARIES})
ENDIF(WIN32)

SET_TARGET_PROPERTIES(LibSL PROPERTIES DEBUG_POSTFIX "-d")
SET_PROPERTY(TARGET LibSL APPEND PROPERTY
   INTERFACE_INCLUDE_DIRECTORIES
 			${CMAKE_CURRENT_SOURCE_DIR}/src/
			${CMAKE_CURRENT_SOURCE_DIR}/src/LibSL
			${CMAKE_CURRENT_SOURCE_DIR}/src/LibSL/loki/
)

ADD_DEFINITIONS(-D_UNICODE -DUNICODE)

# Enable max performance
ADD_DEFINITIONS("-DNDEBUG")
ADD_DEFINITIONS("-DLIBSL_RELEASE")

INSTALL(TARGETS LibSL
  RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
  LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib
  ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/lib
)
