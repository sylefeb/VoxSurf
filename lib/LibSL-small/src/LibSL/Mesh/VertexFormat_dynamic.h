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
#include <LibSL/Errors/Errors.h>
#include <LibSL/Math/Vertex.h>

#include <map>

// -------------------------------------------------
namespace LibSL {
  namespace Mesh {

    // -------------------------------------------------

    template <class MVF_base_attribute> class MVF_get_binding;

    template <class MVF_type> class MVF_get_type;

    // -------------------------------------------------

    class MVF
    {
    public:

      // All bindings
      enum e_Binding {
        Undefined = -1,
        Position  = 0, 
        Normal    = 1, 
        Color0    = 2, 
        Color1    = 3, 
        TexCoord0 = 4, 
        TexCoord1 = 5, 
        TexCoord2 = 6, 
        TexCoord3 = 7, 
        TexCoord4 = 8, 
        TexCoord5 = 9, 
        TexCoord6 = 10, 
        TexCoord7 = 11
      };

      // All types
      enum e_Type {
        Null   = -1,
        Byte   = 0,
        Float  = 1,
        Double = 2,
        Int    = 3
      };

      class Attribute
      {
      public:
        e_Binding binding;       // attribute binding
        e_Type    type;          // type of the attribute
        int       numComponents; // number of components
        int       index;         // index within MVF
        int       offset;        // offset (as bytes) within vertex data field
        int       size_of;       // size of attribute (sizeof(type)*numComponents)

        bool      isValid() const {
          return (binding != Undefined && type != Null && numComponents > -1 
             && index > -1 && offset > -1 && size_of > -1); }
        
        Attribute() : binding(Undefined), type(Null), numComponents(-1), 
                      index(-1), offset(-1), size_of(-1) {}
      };

    protected:

      std::vector<Attribute>   m_Attributes;
      std::map<e_Binding,int>  m_AttributesByBinding;

    public:

      MVF()         { }

      MVF(const MVF *mvf)
      { 
        m_Attributes          = mvf->m_Attributes;
        m_AttributesByBinding = mvf->m_AttributesByBinding;
      }

      template <typename MVF_static> static MVF *make();

      //! Get all attributes (read only)
      const std::vector<Attribute>&   attributes() const { return m_Attributes; }
      //! Finds an attribute by binding (may return NULL if binding is not used)
      const Attribute                *findAttributeByBinding(e_Binding) const;
      //! Adds an attribute (carreful: all attribute fields must be properly filled-in).
      void                            addAttribute(const Attribute& a);
      //! Converts a vertex data into this vertex format
      //    (missing attributes are left untouched)
      void                            convertVertexFrom(void *dst,const void *src,const MVF *srcmvf) const;

      // returns a pointer on a v3f describing the vertex position
      // throws Fatal if position is not present, or not 3 floats
	  // this is the same as attr with the Position binding, just provided for convenience
      LibSL::Math::v3f               *pos3(void *data);
      const LibSL::Math::v3f         *pos3(const void *data) const;

	  // returns a pointer to the start of an attribute - static version
	  // TODO: const version
	  template<typename MVF_attr>
	  typename MVF_attr::type *attr(void *data)
	  {
		  MVF::e_Binding binding = (MVF::e_Binding)MVF_get_binding<typename MVF_attr::base_attribute>::value;
		  MVF::e_Type    type    = (MVF::e_Type)   MVF_get_type   <typename MVF_attr::type          >::value;
		  const Attribute *a = findAttributeByBinding( binding );
		  if (a == NULL) {
		          throw LibSL::Errors::Fatal("MVF::attr - attribute not found in MVF");
		  }
		  sl_assert(a->binding == binding);
		  if (a->numComponents != MVF_attr::components) {
			  throw LibSL::Errors::Fatal("MVF::attr - attribute has incorrect number of components",a->numComponents);
		  }
		  if (a->type != type) {
			  throw LibSL::Errors::Fatal("MVF::attr - attribute type mismatch");
		  }
		  return (typename MVF_attr::type *)((uchar*)data + a->offset);
	  }

	  // returns a pointer to the start of an attribute - dynamic version
	  // TODO: const version
	  void *attr(void *data,MVF::e_Type type,MVF::e_Binding binding,uint nComp);

	  // test whether an attribute is present - static version
	  template<typename MVF_attr>
	  bool hasAttr()
	  {
		  MVF::e_Binding binding = (MVF::e_Binding)MVF_get_binding<typename MVF_attr::base_attribute>::value;
		  MVF::e_Type    type    = (MVF::e_Type)   MVF_get_type   <typename MVF_attr::type          >::value;
		  const Attribute *a = findAttributeByBinding( binding );
		  if (a == NULL) {
			  return false;
		  }
		  sl_assert(a->binding == binding);
		  if (a->numComponents != MVF_attr::components) {
			  return false;
		  }
		  if (a->type != type) {
			  return false;
		  }
		  return true;
	  }

	  // test whether an attribute is present - dynamic version
	  bool hasAttr(MVF::e_Type type,MVF::e_Binding binding,uint nComp);

	  uint sizeOf() const
	  {
		  uint szof = 0;
        ForIndex(a,attributes().size()) {
          szof += attributes()[a].size_of;
        }
        return szof;
      }

      bool load(FILE *f);
      void save(FILE *f) const;

    };

    // -------------------------------------------------
    // Template constructs to convert a static MVF into a dynamic one

    template < > class MVF_get_binding<MVF_BASE_POSITION>  { public: enum { value = MVF::Position};  };
    template < > class MVF_get_binding<MVF_BASE_NORMAL>    { public: enum { value = MVF::Normal};    };
    template < > class MVF_get_binding<MVF_BASE_COLOR0>    { public: enum { value = MVF::Color0};    };
    template < > class MVF_get_binding<MVF_BASE_COLOR1>    { public: enum { value = MVF::Color1};    };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD0> { public: enum { value = MVF::TexCoord0}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD1> { public: enum { value = MVF::TexCoord1}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD2> { public: enum { value = MVF::TexCoord2}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD3> { public: enum { value = MVF::TexCoord3}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD4> { public: enum { value = MVF::TexCoord4}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD5> { public: enum { value = MVF::TexCoord5}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD6> { public: enum { value = MVF::TexCoord6}; };
    template < > class MVF_get_binding<MVF_BASE_TEXCOORD7> { public: enum { value = MVF::TexCoord7}; };

    template < > class MVF_get_type<char>     { public: enum { value = MVF::Byte};   };
    template < > class MVF_get_type<uchar>    { public: enum { value = MVF::Byte};   };
    template < > class MVF_get_type<float>    { public: enum { value = MVF::Float};  };
    template < > class MVF_get_type<double>   { public: enum { value = MVF::Double}; };
    template < > class MVF_get_type<int>      { public: enum { value = MVF::Int};    };
    template < > class MVF_get_type<uint>     { public: enum { value = MVF::Int};    };

    template <class MVF_static,class MVF_current> class MVF_ConvertToDynamic;

    template <class MVF_static> class MVF_ConvertToDynamic<MVF_static,Loki::NullType>
    {
    public:
      MVF_ConvertToDynamic(MVF *) { }
    };

    template <class MVF_static,class Head,class Tail> 
    class MVF_ConvertToDynamic<MVF_static,Loki::Typelist<Head,Tail> > : public MVF_ConvertToDynamic<MVF_static,Tail>
    {
    protected:
      enum { offset  = MVF_offset_item<MVF_static,typename Head::base_attribute>::value };
      enum { index   = MVF_index_item <MVF_static,typename Head::base_attribute>::value };
      enum { size_of = Head::size_of };
    public:
      MVF_ConvertToDynamic(MVF *mvf) : MVF_ConvertToDynamic<MVF_static,Tail>(mvf)
      { 
        MVF::Attribute a;
        a.index         = index;
        a.offset        = offset;
        a.size_of       = size_of;
        a.binding       = (MVF::e_Binding)MVF_get_binding<typename Head::base_attribute>::value;
        a.type          = (MVF::e_Type)   MVF_get_type   <typename Head::type          >::value;
        a.numComponents = Head::components;
        mvf->addAttribute( a );
      }
    };

    /// dynamic from static MVF

    template <typename MVF_static> 
    MVF *MVF::make()
    {
      MVF *mvf = new MVF();
      MVF_ConvertToDynamic<MVF_static,typename Loki::TL::Reverse<MVF_static>::Result > converter(mvf);
      return mvf;
    }

    // -------------------------------------------------
  } // namespace Mesh
} // namespace LibSL
// -------------------------------------------------

std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::e_Binding& b);
std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::e_Type& t);
std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF::Attribute& a);
std::ostream& operator<<(std::ostream& s,const LibSL::Mesh::MVF& mvf);

// -------------------------------------------------
