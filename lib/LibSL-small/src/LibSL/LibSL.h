/* --------------------------------------------------------------------
Author: Sylvain Lefebvre    sylvain.lefebvre@inria.fr

                  Simple Library for Graphics (LibSL)

This software is a computer program whose purpose is to offer a set of
tools to simplify programming real-time computer graphics applications
under OpenGL and DirectX.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.
-------------------------------------------------------------------- */
// ------------------------------------------------------
// LibSL - main include file
// ------------------------------------------------------
//
//
//
// ------------------------------------------------------
// Sylvain Lefebvre - 2006-03-09
// ------------------------------------------------------

#pragma once
#define LIBSL_CORE_INCLUDED

// #pragma message("Including LibSL.h")

#define LIBSL_RELEASE

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif

#include <LibSL/Errors/Errors.h>
#include <LibSL/Memory/Pointer.h>
#include <LibSL/System/System.h>
#include <LibSL/CppHelpers/CppHelpers.h>
#include <LibSL/StlHelpers/StlHelpers.h>
#include <LibSL/SvgHelpers/SvgHelpers.h>

#include <LibSL/Memory/Array.h>
#include <LibSL/Memory/Array2D.h>
#include <LibSL/Memory/Array3D.h>
#include <LibSL/Memory/ArrayTools.h>

#include <LibSL/Math/Math.h>
#include <LibSL/Math/Tuple.h>
#include <LibSL/Math/Vertex.h>
#include <LibSL/Math/Matrix4x4.h>
#include <LibSL/Math/Quaternion.h>

#include <LibSL/Geometry/AAB.h>
#include <LibSL/Mesh/Mesh.h>
#include <LibSL/Mesh/MeshFormat_stl.h>

#include <LibSL/Image/Image.h>
#include <LibSL/Image/Filter.h>
#include <LibSL/Image/ImagePyramid.h>
#include <LibSL/Image/ImageFormat_TGA.h>

// using namespace LibSL;
using namespace LibSL;
using namespace LibSL::Errors;
using namespace LibSL::System;
using namespace LibSL::System::File; 
using namespace LibSL::System::Time;
using namespace LibSL::Memory::Array;
using namespace LibSL::Memory::Pointer;
using namespace LibSL::Math; 
using namespace LibSL::Geometry;
using namespace LibSL::CppHelpers;
using namespace LibSL::StlHelpers;
using namespace LibSL::Mesh;
using namespace LibSL::Image;
using namespace LibSL::Filter;

#define LIBSL_WIN32_FIX                \
LibSL::Image::ImageFormat_TGA   s_TGA; \
LibSL::Mesh::MeshFormat_stl     s_Stl;

#define LIBSL_STATIC_FIX LIBSL_WIN32_FIX

// ------------------------------------------------------
