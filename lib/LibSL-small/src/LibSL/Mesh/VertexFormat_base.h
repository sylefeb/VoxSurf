/* --------------------------------------------------------------------
Author: Sylvain Lefebvre    sylvain.lefebvre@sophia.inria.fr

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

#pragma once

namespace LibSL {
  namespace Mesh {

// Base attributes
class MVF_BASE_POSITION  {};
class MVF_BASE_NORMAL    {};
class MVF_BASE_COLOR0    {};
class MVF_BASE_COLOR1    {};
class MVF_BASE_TEXCOORD0 {};
class MVF_BASE_TEXCOORD1 {};
class MVF_BASE_TEXCOORD2 {};
class MVF_BASE_TEXCOORD3 {};
class MVF_BASE_TEXCOORD4 {};
class MVF_BASE_TEXCOORD5 {};
class MVF_BASE_TEXCOORD6 {};
class MVF_BASE_TEXCOORD7 {};
class MVF_BASE_BINORMAL     {}; /// TODO/FIXME: these are not sent as vertex attributes ... remove?
class MVF_BASE_TANGENT      {};
class MVF_BASE_BONE_INDICES {};
class MVF_BASE_BONE_WEIGHTS {};

// Attribute template
template <class Base,class StorageType,unsigned int sz> 
class MVF_ITEM
{
public:
  typedef Base        base_attribute;
  typedef StorageType type;
  enum               {size_of = sizeof(StorageType)*sz};
  enum               {components = sz};

  StorageType         data[sz];
};

// Empty attribute
typedef char MVF_null;
class   MVF_empty : public MVF_ITEM<Loki::NullType,MVF_null,1> {};

// Attributes
template <class StorageType,unsigned int sz> 
  class MVF_POSITION : public MVF_ITEM<MVF_BASE_POSITION,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_NORMAL : public MVF_ITEM<MVF_BASE_NORMAL,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_COLOR0 : public MVF_ITEM<MVF_BASE_COLOR0,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_COLOR1 : public MVF_ITEM<MVF_BASE_COLOR1,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD0 : public MVF_ITEM<MVF_BASE_TEXCOORD0,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD1 : public MVF_ITEM<MVF_BASE_TEXCOORD1,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD2 : public MVF_ITEM<MVF_BASE_TEXCOORD2,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD3 : public MVF_ITEM<MVF_BASE_TEXCOORD3,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD4 : public MVF_ITEM<MVF_BASE_TEXCOORD4,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD5 : public MVF_ITEM<MVF_BASE_TEXCOORD5,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD6 : public MVF_ITEM<MVF_BASE_TEXCOORD6,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TEXCOORD7 : public MVF_ITEM<MVF_BASE_TEXCOORD7,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_BINORMAL : public MVF_ITEM<MVF_BASE_BINORMAL,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_TANGENT : public MVF_ITEM<MVF_BASE_TANGENT,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_BONE_INDICES : public MVF_ITEM<MVF_BASE_BONE_INDICES,StorageType,sz> {};
template <class StorageType,unsigned int sz> 
  class MVF_BONE_WEIGHTS : public MVF_ITEM<MVF_BASE_BONE_WEIGHTS,StorageType,sz> {};


// -------------------------------------------------
// typedef some standard attributes
// -------------------------------------------------

// -> These correspond directly to GPU bindings (ie. they are understood by GPUMesh)
typedef MVF_POSITION<float,2>       mvf_vertex_2f;
typedef MVF_POSITION<double,2>      mvf_vertex_2d;
typedef MVF_POSITION<float,3>       mvf_vertex_3f;
typedef MVF_POSITION<double,3>      mvf_vertex_3d;
typedef MVF_POSITION<float,4>       mvf_vertex_4f;
typedef MVF_POSITION<double,4>      mvf_vertex_4d;
typedef MVF_POSITION<float,2>       mvf_position_2f;
typedef MVF_POSITION<double,2>      mvf_position_2d;
typedef MVF_POSITION<float,3>       mvf_position_3f;
typedef MVF_POSITION<double,3>      mvf_position_3d;
typedef MVF_POSITION<float,4>       mvf_position_4f;
typedef MVF_POSITION<double,4>      mvf_position_4d;
typedef MVF_NORMAL<float,3>         mvf_normal_3f;
typedef MVF_NORMAL<double,3>        mvf_normal_3d;
typedef MVF_COLOR0<float,3>         mvf_color0_3f;
typedef MVF_COLOR0<double,3>        mvf_color0_3d;
typedef MVF_COLOR0<float,4>         mvf_color0_4f;
typedef MVF_COLOR0<double,4>        mvf_color0_4d;
typedef MVF_COLOR0<unsigned char,4> mvf_color0_rgba;
typedef MVF_COLOR0<unsigned char,3> mvf_color0_rgb;
typedef MVF_COLOR1<float,3>         mvf_color1_3f;
typedef MVF_COLOR1<double,3>        mvf_color1_3d;
typedef MVF_COLOR1<float,4>         mvf_color1_4f;
typedef MVF_COLOR1<double,4>        mvf_color1_4d;
typedef MVF_COLOR1<unsigned char,4> mvf_color1_rgba;
typedef MVF_COLOR1<unsigned char,3> mvf_color1_rgb;
typedef MVF_TEXCOORD0<float,1>      mvf_texcoord0_1f;
typedef MVF_TEXCOORD0<double,1>     mvf_texcoord0_1d;
typedef MVF_TEXCOORD0<float,2>      mvf_texcoord0_2f;
typedef MVF_TEXCOORD0<double,2>     mvf_texcoord0_2d;
typedef MVF_TEXCOORD0<float,3>      mvf_texcoord0_3f;
typedef MVF_TEXCOORD0<double,3>     mvf_texcoord0_3d;
typedef MVF_TEXCOORD0<float,4>      mvf_texcoord0_4f;
typedef MVF_TEXCOORD0<double,4>     mvf_texcoord0_4d;
typedef MVF_TEXCOORD1<float,1>      mvf_texcoord1_1f;
typedef MVF_TEXCOORD1<double,1>     mvf_texcoord1_1d;
typedef MVF_TEXCOORD1<float,2>      mvf_texcoord1_2f;
typedef MVF_TEXCOORD1<double,2>     mvf_texcoord1_2d;
typedef MVF_TEXCOORD1<float,3>      mvf_texcoord1_3f;
typedef MVF_TEXCOORD1<double,3>     mvf_texcoord1_3d;
typedef MVF_TEXCOORD1<float,4>      mvf_texcoord1_4f;
typedef MVF_TEXCOORD1<double,4>     mvf_texcoord1_4d;
typedef MVF_TEXCOORD2<float,1>      mvf_texcoord2_1f;
typedef MVF_TEXCOORD2<double,1>     mvf_texcoord2_1d;
typedef MVF_TEXCOORD2<float,2>      mvf_texcoord2_2f;
typedef MVF_TEXCOORD2<double,2>     mvf_texcoord2_2d;
typedef MVF_TEXCOORD2<float,3>      mvf_texcoord2_3f;
typedef MVF_TEXCOORD2<double,3>     mvf_texcoord2_3d;
typedef MVF_TEXCOORD2<float,4>      mvf_texcoord2_4f;
typedef MVF_TEXCOORD2<double,4>     mvf_texcoord2_4d;
typedef MVF_TEXCOORD3<float,1>      mvf_texcoord3_1f;
typedef MVF_TEXCOORD3<double,1>     mvf_texcoord3_1d;
typedef MVF_TEXCOORD3<float,2>      mvf_texcoord3_2f;
typedef MVF_TEXCOORD3<double,2>     mvf_texcoord3_2d;
typedef MVF_TEXCOORD3<float,3>      mvf_texcoord3_3f;
typedef MVF_TEXCOORD3<double,3>     mvf_texcoord3_3d;
typedef MVF_TEXCOORD3<float,4>      mvf_texcoord3_4f;
typedef MVF_TEXCOORD3<double,4>     mvf_texcoord3_4d;
typedef MVF_TEXCOORD4<float,1>      mvf_texcoord4_1f;
typedef MVF_TEXCOORD4<double,1>     mvf_texcoord4_1d;
typedef MVF_TEXCOORD4<float,2>      mvf_texcoord4_2f;
typedef MVF_TEXCOORD4<double,2>     mvf_texcoord4_2d;
typedef MVF_TEXCOORD4<float,3>      mvf_texcoord4_3f;
typedef MVF_TEXCOORD4<double,3>     mvf_texcoord4_3d;
typedef MVF_TEXCOORD4<float,4>      mvf_texcoord4_4f;
typedef MVF_TEXCOORD4<double,4>     mvf_texcoord4_4d;
typedef MVF_TEXCOORD5<float,1>      mvf_texcoord5_1f;
typedef MVF_TEXCOORD5<double,1>     mvf_texcoord5_1d;
typedef MVF_TEXCOORD5<float,2>      mvf_texcoord5_2f;
typedef MVF_TEXCOORD5<double,2>     mvf_texcoord5_2d;
typedef MVF_TEXCOORD5<float,3>      mvf_texcoord5_3f;
typedef MVF_TEXCOORD5<double,3>     mvf_texcoord5_3d;
typedef MVF_TEXCOORD5<float,4>      mvf_texcoord5_4f;
typedef MVF_TEXCOORD5<double,4>     mvf_texcoord5_4d;
typedef MVF_TEXCOORD6<float,1>      mvf_texcoord6_1f;
typedef MVF_TEXCOORD6<double,1>     mvf_texcoord6_1d;
typedef MVF_TEXCOORD6<float,2>      mvf_texcoord6_2f;
typedef MVF_TEXCOORD6<double,2>     mvf_texcoord6_2d;
typedef MVF_TEXCOORD6<float,3>      mvf_texcoord6_3f;
typedef MVF_TEXCOORD6<double,3>     mvf_texcoord6_3d;
typedef MVF_TEXCOORD6<float,4>      mvf_texcoord6_4f;
typedef MVF_TEXCOORD6<double,4>     mvf_texcoord6_4d;
typedef MVF_TEXCOORD7<float,1>      mvf_texcoord7_1f;
typedef MVF_TEXCOORD7<double,1>     mvf_texcoord7_1d;
typedef MVF_TEXCOORD7<float,2>      mvf_texcoord7_2f;
typedef MVF_TEXCOORD7<double,2>     mvf_texcoord7_2d;
typedef MVF_TEXCOORD7<float,3>      mvf_texcoord7_3f;
typedef MVF_TEXCOORD7<double,3>     mvf_texcoord7_3d;
typedef MVF_TEXCOORD7<float,4>      mvf_texcoord7_4f;
typedef MVF_TEXCOORD7<double,4>     mvf_texcoord7_4d;

// -> These are defined for convenience but do not correspond to GPU bindings (GPUMesh will reject them)
typedef MVF_BINORMAL<float,3>       mvf_binormal_3f;
typedef MVF_BINORMAL<double,3>      mvf_binormal_3d;
typedef MVF_TANGENT<float,3>        mvf_tangent_3f;
typedef MVF_TANGENT<double,3>       mvf_tangent_3d;
typedef MVF_BONE_INDICES<int,1>     mvf_bone_indices_1i;
typedef MVF_BONE_INDICES<int,2>     mvf_bone_indices_2i;
typedef MVF_BONE_INDICES<int,3>     mvf_bone_indices_3i;
typedef MVF_BONE_INDICES<int,4>     mvf_bone_indices_4i;
typedef MVF_BONE_WEIGHTS<double,1>  mvf_bone_weights_1d;
typedef MVF_BONE_WEIGHTS<double,2>  mvf_bone_weights_2d;
typedef MVF_BONE_WEIGHTS<double,3>  mvf_bone_weights_3d;
typedef MVF_BONE_WEIGHTS<double,4>  mvf_bone_weights_4d;
typedef MVF_BONE_WEIGHTS<float,1>   mvf_bone_weights_1f;
typedef MVF_BONE_WEIGHTS<float,2>   mvf_bone_weights_2f;
typedef MVF_BONE_WEIGHTS<float,3>   mvf_bone_weights_3f;
typedef MVF_BONE_WEIGHTS<float,4>   mvf_bone_weights_4f;

    // -------------------------------------------------
  } // namespace Mesh
} // namespace LibSL
// -------------------------------------------------
