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
// ----------------------------------------------------------------------------
// VertexFormat.h
//
//    This file defines the tools needed to specify and work with 
//    API independent vertex format. The vertex format is defined
//    as a typelist. Each type of this list is defining an attribute 
//    (ex: the second texture coordinate as 2 floats is the 
//     mvf_texcoord1_2f type - or MVF_TEXCOORD1<float,2>).
//    Each attribute is associated to a base attribute. 
//    (ex: mvf_texcoord1_2f and mvf_texcoord1_3d are both associated
//     with MVF_BASE_TEXCOORD1 because they both correspond to a definition
//    of the second texture coordinate).
//    Given an attribute T, the base attribute is T::base_attribute
//
//    A format can easily be described by using the FormatDescriptor
//    class templated by the vertex format. It allows to know which
//    base attributes are present, the size of each attribute, the 
//    size of a vertex with this format, the offset of an attribute
//    in a table of bytes representing a vertex.
//
//
//                                                    (c) Sylvain Lefebvre 2003
// ----------------------------------------------------------------------------

#pragma once

#include <loki/Typelist.h>
#include <LibSL/System/Types.h>

// -------------------------------------------------
// Vertex attributes
// -------------------------------------------------
#include "VertexFormat_base.h"

// -------------------------------------------------
namespace LibSL {
  namespace Mesh {
    // -------------------------------------------------


    // -------------------------------------------------
    // Sizeof vertex attributes and vertex format
    // -------------------------------------------------

    template <class TList> class MVF_sizeof;

    template <> class MVF_sizeof< Loki::NullType > {public: enum {value = 0}; };

    template <class Head,class Tail> class MVF_sizeof< Loki::Typelist<Head,Tail> >
    { 
    private:

    public:

      enum {value = Head::size_of + MVF_sizeof<Tail>::value};
    };

    // -------------------------------------------------
    // Check whether an attribute corresponds to an item
    // -------------------------------------------------

    template <class T,class Item> class MVF_is_item
    {
    public:
      enum {value = 0};
    };

    template <class T> class MVF_is_item<T,typename T::base_attribute>
    {
    public:
      enum {value = 1};
    };

    // -------------------------------------------------
    // Vertex attribute offset
    // -------------------------------------------------

    template <class TList,class T> class MVF_offset;

    template <class T> class MVF_offset<Loki::NullType,T>
    {
    public:
      enum { value = -1 }; 
    };

    template <class Tail,class T> class MVF_offset<Loki::Typelist<T,Tail>,T>
    {
    public:
      enum { value = 0 };   
    };

    template <class Head,class Tail,class T> class MVF_offset<Loki::Typelist<Head,Tail>,T>
    {
    public:
      enum { temp  = (MVF_offset<Tail,T>::value) };
      enum { value = ((temp == -1) ? -1 : (Head::size_of + temp))};
    };

    // -------------------------------------------------
    // Vertex attribute offset from item
    // -------------------------------------------------

    template <class TList,class Item> class MVF_offset_item;

    template <class Item> class MVF_offset_item<Loki::NullType,Item>
    {
    public:
      enum {value = -1};
    };

    template <class Head,class Tail,class Item> class MVF_offset_item<Loki::Typelist<Head,Tail>,Item>
    {
    private:
      enum { stop = MVF_is_item<Head,Item>::value     };
      enum { next = MVF_offset_item<Tail,Item>::value };
    public:
      enum { value = (stop == 0) ? ((next == -1) ? -1 : (Head::size_of + next)) : 0 };
    };

    // -------------------------------------------------
    // Retrieve index of an item
    // -------------------------------------------------

    template <class TList,class Item> class MVF_index_item;

    template <class Item> class MVF_index_item<Loki::NullType,Item>
    {
    public:
      enum {value = -1};
    };

    template <class Head,class Tail,class Item> class MVF_index_item<Loki::Typelist<Head,Tail>,Item>
    {
    private:
      enum { stop = MVF_is_item<Head,Item>::value    };
      enum { next = MVF_index_item<Tail, Item>::value};
    public:
      enum { value = (stop == 0) ? ((next == -1) ? -1 : (1+next)) : 0 };
    };

    // -------------------------------------------------
    // CompatibleFormat
    //   test if MVF1 is compatible with MVF2
    //   this is the case if MVF2 is a superset of MVF1
    // -------------------------------------------------

    template <class MVF1,class MVF2> class CompatibleFormat;

    template <class MVF2> class CompatibleFormat<Loki::NullType,MVF2>
    {
    public:
      enum {value = 1};
    };

    template <class Head,class Tail,class MVF2> class CompatibleFormat<Loki::Typelist<Head,Tail>,MVF2>
    {
    private:
      enum { temp = (CompatibleFormat<Tail,MVF2>::value) };
    public:
      enum { value = temp * ((MVF_index_item<MVF2,typename Head::base_attribute>::value > -1)? 1 : 0)};
    };

    // -------------------------------------------------
    //  define a macro helping to check format compatibility

#define CHECK_COMPATIBLE_FORMAT(mvf1,mvf2) { \
  const int __d=CompatibleFormat<mvf1,mvf2>::value; \
  LOKI_STATIC_CHECK(__d,vertex_format_not_compatible); }

    // -------------------------------------------------
    // Converts data corresponding to a format into another, compatible one
    // -------------------------------------------------

    template <class MVF_Dst,class MVF_Src,class MVF_Current> class ConvertToFormat_rec;

    template <class MVF_Dst,class MVF_Src> class ConvertToFormat_rec<MVF_Dst,MVF_Src,Loki::NullType>
    {
    public:
      ConvertToFormat_rec(uchar *dst,const uchar *src) { }
    };

    template <class MVF_Dst,class MVF_Src,class Head,class Tail> 
    class ConvertToFormat_rec<MVF_Dst,MVF_Src,Loki::Typelist<Head,Tail> > : public ConvertToFormat_rec<MVF_Dst,MVF_Src,Tail>
    {
    private:
      enum { offset_dst = MVF_offset_item<MVF_Dst,typename Head::base_attribute>::value };
      enum { offset_src = MVF_offset_item<MVF_Src,typename Head::base_attribute>::value };
      enum { size_of    = Head::size_of };
    public:
      ConvertToFormat_rec(uchar *dst,const uchar *src) : ConvertToFormat_rec<MVF_Dst,MVF_Src,Tail>(dst,src)
      { 
        memcpy(dst + offset_dst , src + offset_src , size_of);
      }
    };

    template <class MVF_Dst,class MVF_Src> 
    class ConvertToFormat : public ConvertToFormat_rec<MVF_Dst,MVF_Src,MVF_Src>
    {
    public:
      ConvertToFormat(void *dst,const void *src) : ConvertToFormat_rec<MVF_Dst,MVF_Src,MVF_Src>((uchar*)dst,(const uchar*)src)
      { 

      }
    };

    // -------------------------------------------------
  } // namespace Mesh
} // namespace LibSL
// -------------------------------------------------
