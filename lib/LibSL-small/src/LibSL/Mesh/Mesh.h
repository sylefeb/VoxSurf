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
// ------------------------------------------------------
// LibSL::Mesh
// ------------------------------------------------------
//
// Mesh generic class
//
// T_VertexData parameter format:
//  - required to have a member 'pos' of type v3f
//
// ------------------------------------------------------
// Sylvain Lefebvre - 2006-11-15
// ------------------------------------------------------

#pragma once

#include <LibSL/LibSL.common.h>
#include <LibSL/Errors/Errors.h>
#include <LibSL/Math/Matrix4x4.h>
#include <LibSL/Math/Vertex.h>
#include <LibSL/Memory/Pointer.h>
#include <LibSL/Memory/Array.h>
#include <LibSL/Math/Vertex.h>
#include <LibSL/Geometry/AAB.h>
#include <LibSL/Mesh/VertexFormat.h>
#include <loki/static_check.h>

#include <cstring>

// ------------------------------------------------------

namespace LibSL {
  namespace Mesh {

    /// Forward declaration of generic class
    template <class T_VertexData>    class TriangleMesh_generic;

    /// Forward declaration of vertex container class
    template <typename T_VertexData> class VertexContainer;

    /// Base Triangle Mesh class
    class TriangleMesh
      //        : public LibSL::Memory::TraceLeaks::LeakProbe<TriangleMesh>
    {
    public:

      typedef LibSL::Math::v3f           v3;
      typedef LibSL::Math::v4f           v4;
      typedef LibSL::Math::Tuple<uint,3> t_Triangle;

      // struct to hold surface information
      typedef struct {
        std::string                       textureName;  // texture name for the surface
        v3                                diffuse;      // diffuse color for the surface
        LibSL::Memory::Array::Array<uint> triangleIds;  // ids of triangles using the texture
      } t_SurfaceNfo;

    protected:

      virtual void truncateVertices(uint sz)                                                  =0;
      virtual void reorderVerticesAndTruncate(const LibSL::Memory::Array::Array<uint>& order) = 0;

      bool                   m_BBoxComputed;
      LibSL::Geometry::AABox m_BBox;

    public:

      typedef LibSL::Memory::Pointer::AutoPtr<TriangleMesh> t_AutoPtr;

      TriangleMesh() : m_BBoxComputed(false) {}
      virtual ~TriangleMesh() {}

      virtual void                allocate(uint numvert,uint numtris,uint numsurfaces=0) =0;

      virtual uint                numVertices()         const   =0;
      virtual const v3&           posAt(uint n)         const   =0;
      virtual v3&                 posAt(uint n)                 =0;
      virtual const void         *vertexDataAt(uint n)  const   =0;
      virtual void               *vertexDataAt(uint n)          =0;
      virtual uint                sizeOfVertexData()    const   =0;

      virtual uint                numTriangles()        const   =0;
      virtual const t_Triangle&   triangleAt(uint t)    const   =0;
      virtual t_Triangle&         triangleAt(uint t)            =0;

      virtual const LibSL::Memory::Pointer::AutoPtr<MVF>& mvf()                 const   =0;
      virtual void                setMvf(const LibSL::Memory::Pointer::AutoPtr<MVF>&)   =0;

      //! information about textured surfaces
      virtual uint                 numSurfaces()                      const =0;
      virtual const t_SurfaceNfo&  surfaceAt(uint s)                  const =0;
      virtual       t_SurfaceNfo&  surfaceAt(uint s)                        =0;
      virtual std::string&         surfaceTextureName (uint s)              =0;
      virtual const std::string&   surfaceTextureName (uint s)        const =0;
      virtual uint                 surfaceNumTriangles(uint s)        const =0;
      virtual uint                 surfaceTriangleIdAt(uint s,uint t) const =0;

      //! return bounding box
      const LibSL::Geometry::AABox& bbox();
      LibSL::Geometry::AABox computeBBox();

      //! scale mesh into unit cube
      void scaleToUnitCube(float scale=1.0f,bool center=true);

      //! apply transform to mesh
      void applyTransform(const Math::m4x4f& trsf);

      //! center mesh on a specific point
      LibSL::Math::v3f  centerOn(const LibSL::Math::v3f& ctr);

      //! merge vertices that are close to each other
      //  vertices closer than radius will be snapped together
      //  (gridsz controls the resolution of the acceleration grid)
      void mergeVertices(float radius=1e-9f,uint gridsz=256);
      void mergeVerticesExact();

      //! create a new mesh wich is the same but with no shared vertices between triangles
      virtual TriangleMesh *splitTriangles() const;

      //! reorient triangles consistently
      //  only succeeds if two-manifold
      virtual void reorientTriangles();

      //! swap axes (specify a 3 character chain for reordering "yxz" ... A capital letter negates the component)
      void swapAxes(const char *swizzle="xyz");

      //! instanciate a mesh using specialized class
      virtual TriangleMesh *newInstance()  const = 0;

      //! clone the current mesh
      virtual TriangleMesh *clone()        const;

      //! dynamic cast - creates a new mesh with a different vertex format
      //               - only necessary data is copied, missing entries are zeroed
      //               - contrary to static cast, the conversion is performed using dynamic
      //                 vertex formats given as parameters
      //               - attributes may not match exactly: mvf_texcoord0_2f will get copied
      //                 into mvf_texcoord0_3f (this is different from staticCast)
      //               - NOTE: a staticCast is available in TriangleMesh_generic
      //  T_VertexDataDst:   destination data structure
      //  mvfdst         :   destination vertex format
      //  mvfsrc         :   source vertex format (how to interpret this mesh vertex format)
      template <class T_VertexDataDst>
      TriangleMesh_generic<T_VertexDataDst> *dynamicCast(const MVF *mvfdst,const MVF *mvfsrc) const
      {
        sl_assert(mvfdst != NULL);
        sl_assert(mvfsrc != NULL);
        if (mvfsrc->sizeOf() != sizeOfVertexData()) {
          throw Errors::Fatal("Source vertex format does not have mesh vertex data_size\n");
        }
        if (mvfdst->sizeOf() != VertexContainer<T_VertexDataDst>::size_of && VertexContainer<T_VertexDataDst>::size_of > -1) {
          throw Errors::Fatal("Destination vertex format does not have mesh vertex data_size\n");
        }
        TriangleMesh_generic<T_VertexDataDst> *mesh = new TriangleMesh_generic<T_VertexDataDst>(LibSL::Memory::Pointer::AutoPtr<MVF>(new MVF(mvfdst)));
        // allocate
        mesh->allocate(numVertices(),numTriangles(),numSurfaces());
        // copy triangles and surfaces
        ForIndex(t,numTriangles()) {
          mesh->triangleAt(t) = triangleAt(t);
        }
        ForIndex(s,numSurfaces()) {
          mesh->surfaceAt(s)  = surfaceAt(s);
        }
        // fill-in with data, converting formats
        ForIndex(v,numVertices()) {
          memset(mesh->vertexDataAt(v),0,mesh->sizeOfVertexData());
          mvfdst->convertVertexFrom( mesh->vertexDataAt(v), vertexDataAt(v),mvfsrc );
        }
        return mesh;
      }

      // All instanciations of generic class are friends
      template <class T> friend class TriangleMesh_generic;
    };

    /// Helper pointer onto base class
    typedef LibSL::Memory::Pointer::AutoPtr<TriangleMesh> TriangleMesh_Ptr;

    /// Vertex containers
    //  Default
    template <typename T_VertexData>
    class VertexContainer : public T_VertexData
    {
    public:
      enum { requiresMVF  = 0 };
      enum { size_of      = sizeof(T_VertexData) };
      void                     init(MVF *mvf)         {  }
      T_VertexData            *getRawData()           { return  static_cast<      T_VertexData *>(this); }
      const T_VertexData      *getRawData()     const { return  static_cast<const T_VertexData *>(this); }
      const TriangleMesh::v3&  getPos(MVF *mvf) const { return  this->pos; }
      TriangleMesh::v3&        getPos(MVF *mvf)       { return  this->pos; }

      VertexContainer()                          { }
      VertexContainer(const T_VertexData& v)     { *static_cast<T_VertexData *>(this) = v; }
    };
    // LibSL::Array
    template <typename T_Type> class VertexContainer < LibSL::Memory::Array::Array<T_Type> > : public LibSL::Memory::Array::Array<T_Type>
    {
    private:
      typedef LibSL::Memory::Array::Array<T_Type> t_Base;
    public:
      enum { requiresMVF  =  1 };
      enum { size_of      = -1 }; // automatically adapts to MVF
      void                     init(MVF *mvf)         { t_Base::allocate( mvf->sizeOf() / sizeof(T_Type) ); }
      T_Type                  *getRawData()           { sl_assert(!t_Base::empty()); return t_Base::raw(); }
      const T_Type            *getRawData()     const { sl_assert(!t_Base::empty()); return t_Base::raw(); }
      const TriangleMesh::v3&  getPos(MVF *mvf) const { sl_assert(!t_Base::empty()); return (*mvf->pos3( t_Base::raw() )); }
      TriangleMesh::v3&        getPos(MVF *mvf)       { sl_assert(!t_Base::empty()); return (*mvf->pos3( t_Base::raw() )); }
    };
    // std::vector
    template <typename T_Type> class VertexContainer < std::vector<T_Type> > : public std::vector<T_Type>
    {
    private:
      typedef std::vector<T_Type> t_Base;
    public:
      enum { requiresMVF  =  1 };
      enum { size_of      = -1 }; // automatically adapts to MVF
      void                     init(MVF *mvf)         { t_Base::resize( mvf->sizeOf() / sizeof(T_Type) ); }
      T_Type                  *getRawData()           { sl_assert(!t_Base::empty()); return &((*this)[0]); }
      const T_Type            *getRawData()     const { sl_assert(!t_Base::empty()); return &((*this)[0]); }
      const TriangleMesh::v3&  getPos(MVF *mvf) const { sl_assert(!t_Base::empty()); return (*mvf->pos3( &((*this)[0]) )); }
      TriangleMesh::v3&        getPos(MVF *mvf)       { sl_assert(!t_Base::empty()); return (*mvf->pos3( &((*this)[0]) )); }
    };

    /// Generic TriangleMesh class
    template <class T_VertexData>
    class TriangleMesh_generic : public TriangleMesh
    {
    public:

      typedef VertexContainer<T_VertexData> t_Vertex;

    protected:

      LibSL::Memory::Array::Array<t_Vertex>      m_Vertices;
      LibSL::Memory::Array::Array<t_Triangle>    m_Triangles;
      LibSL::Memory::Array::Array<t_SurfaceNfo>  m_Surfaces;
      LibSL::Geometry::AAB<3>                    m_Box;
      LibSL::Memory::Pointer::AutoPtr<MVF>       m_MVF;

    protected:

      void allocateVertices(uint numvert)
      {
		  sl_assert( numvert > 0 );
        m_Vertices  .allocate(numvert);
        if (t_Vertex::requiresMVF) {
          ForArray(m_Vertices,v) {
            m_Vertices[v].init(m_MVF.raw());
          }
        }
      }

      virtual void truncateVertices(uint sz)
      {
        if (sz == m_Vertices.size()) {
          return;
        }
        m_Vertices.truncate(sz);
      }

      virtual void reorderVerticesAndTruncate(const LibSL::Memory::Array::Array<uint>& order)
      {
        LibSL::Memory::Array::Array<t_Vertex> vertices(order.size());
        ForArray(order,o) {
          vertices[o] = m_Vertices[order[o]];
        }
        m_Vertices = vertices;
      }

    public:

      TriangleMesh_generic(LibSL::Memory::Pointer::AutoPtr<MVF> mvf = LibSL::Memory::Pointer::AutoPtr<MVF>())
      {
        if (t_Vertex::requiresMVF) { sl_assert(!mvf.isNull()); }
        setMvf(mvf);
      }

      TriangleMesh_generic(uint numvert, uint numtris, uint numsurfaces = 0, LibSL::Memory::Pointer::AutoPtr<MVF> mvf = LibSL::Memory::Pointer::AutoPtr<MVF>())
      {
        if (t_Vertex::requiresMVF) { sl_assert(!mvf.isNull()); }
        setMvf(mvf);
        // allocate data
        allocate(numvert,numtris,numsurfaces);
      }

      virtual void allocate(uint numvert,uint numtris,uint numsurfaces)
      {
        allocateVertices(numvert);
		sl_assert( numtris > 0 );
        m_Triangles .allocate(numtris);
        if (numsurfaces > 0) {
          m_Surfaces.allocate(numsurfaces);
        }
      }

      const T_VertexData&        vertexAt(uint n)    const
      {
        return (m_Vertices[n]);
      }

      T_VertexData&              vertexAt(uint n)
      { return (m_Vertices[n]); }

      virtual const v3&          posAt(uint n)       const
      { return (m_Vertices[n].getPos(mvf().raw())); }

      virtual v3&                posAt(uint n)
      { return (m_Vertices[n].getPos(mvf().raw())); }

      virtual const t_Triangle&  triangleAt(uint t)  const
      { return (m_Triangles[t]); }

      virtual       t_Triangle&  triangleAt(uint t)
      { return (m_Triangles[t]); }

      virtual const t_SurfaceNfo&  surfaceAt(uint s)  const
      { return (m_Surfaces[s]); }

      virtual       t_SurfaceNfo&  surfaceAt(uint s)
      { return (m_Surfaces[s]); }

      virtual uint              numVertices()        const
      { return (m_Vertices.size());  }

      virtual uint              numTriangles()       const
      { return (m_Triangles.size()); }

      const LibSL::Memory::Array::Array<t_Vertex>&     vertices()  const
      { return (m_Vertices); }

      const LibSL::Memory::Array::Array<t_Triangle>&   triangles() const
      { return (m_Triangles); }

      const LibSL::Memory::Array::Array<t_SurfaceNfo>& surfaces()  const
      { return (m_Surfaces); }

      virtual const void       *vertexDataAt(uint n) const
      { return (const void *)( m_Vertices[n].getRawData() ); }

      virtual       void       *vertexDataAt(uint n)
      { return (      void *)( m_Vertices[n].getRawData() ); }

      virtual uint              sizeOfVertexData()   const
      {
        if (t_Vertex::requiresMVF) {
          return mvf()->sizeOf();
        } else {
          return (t_Vertex::size_of);
        }
      }

      virtual uint              numSurfaces()                      const
      {
        return (m_Surfaces.size());
      }

      virtual std::string&      surfaceTextureName (uint s)
      {
        return (m_Surfaces[s].textureName);
      }

      virtual const std::string& surfaceTextureName (uint s)    const
      {
        return (m_Surfaces[s].textureName);
      }

      virtual uint              surfaceNumTriangles(uint s)        const
      {
        return (m_Surfaces[s].triangleIds.size());
      }

      virtual uint              surfaceTriangleIdAt(uint s,uint t) const
      {
        return (m_Surfaces[s].triangleIds[t]);
      }

      virtual const LibSL::Memory::Pointer::AutoPtr<MVF>& mvf() const
      {
        return m_MVF;
      }

      virtual void                setMvf(const LibSL::Memory::Pointer::AutoPtr<MVF>& mvf)
      {
        if (mvf.isNull()) {
          m_MVF = LibSL::Memory::Pointer::AutoPtr<MVF>();
          return;
        }
        m_MVF = mvf;
        // sanity check
        if (mvf->sizeOf() != sizeOfVertexData()) {
          throw Errors::Fatal("TriangleMesh_generic::setMvf - size of MVF (%d) must match size of vertex data (%d)",mvf->sizeOf(),sizeOfVertexData());
        }
      }

      virtual TriangleMesh *newInstance()            const
      { return (new TriangleMesh_generic<T_VertexData>()); }


      //! static cast - creates a new mesh with a different vertex format
      //              - only necessary data is copied, missing entries are zeroed
      //              - everything is resolved at compile time
      //              - attributes must match exactly: mvf_texcoord0_2f will not get copied
      //                into mvf_texcoord0_3f (see dynamicCast if you need this)
      //  T_VertexDataDst:   destination data structure
      //  T_VertexFormatDst: destination format description
      //  T_VertexFormatSrc: source format description (which sizeof must be equal to this mesh data sizeof)
      template <class T_VertexDataDst,class T_VertexFormatDst,class T_VertexFormatSrc>
      TriangleMesh_generic<T_VertexDataDst> *staticCast() const
      {
        LOKI_STATIC_CHECK((int)MVF_sizeof<T_VertexFormatSrc>::value == VertexContainer<T_VertexData   >::size_of || t_Vertex::size_of == -1 ,source_vertex_format_does_not_have_mesh_vertex_data_size);
        LOKI_STATIC_CHECK((int)MVF_sizeof<T_VertexFormatDst>::value == VertexContainer<T_VertexDataDst>::size_of || t_Vertex::size_of == -1 ,destination_vertex_format_does_not_have_mesh_vertex_data_size);
        TriangleMesh_generic<T_VertexDataDst> *mesh( new TriangleMesh_generic<T_VertexDataDst>( LibSL::Memory::Pointer::AutoPtr<MVF>(MVF::make<T_VertexFormatDst>()) ) );
        // copy triangles and surfaces
        mesh->m_Triangles = m_Triangles;
        mesh->m_Surfaces  = m_Surfaces;
        // allocate vertices
        mesh->allocateVertices(numVertices());
        // fill-in with data, converting formats
        ForIndex(v,numVertices()) {
          memset(mesh->vertexDataAt(v),0,mesh->sizeOfVertexData());
          ConvertToFormat<T_VertexFormatDst,T_VertexFormatSrc> converter( mesh->vertexDataAt(v) , vertexDataAt(v) );
        }
        return mesh;
      }

      // All instanciations are friends
      template <class T> friend class TriangleMesh_generic;
    };


    //// Loading and saving meshes

    /// TriangleMesh format plugin (abstract)
    class TriangleMeshFormat_plugin
      //  : public LibSL::Memory::TraceLeaks::LeakProbe<TriangleMeshFormat_plugin>
    {
    public:
      virtual void          save(const char *,const TriangleMesh *) const =0;
      virtual TriangleMesh *load(const char *)                      const =0;
      virtual const char   *signature()                             const =0;
      virtual ~TriangleMeshFormat_plugin() {}
    };

    /// TriangleMesh format manager (singleton)
    class TriangleMeshFormatManager
    {
    private:

      std::map<std::string,TriangleMeshFormat_plugin*> m_Plugins;

      static TriangleMeshFormatManager *s_Manager;
      TriangleMeshFormatManager();

    public:

      ~TriangleMeshFormatManager();

      void registerPlugin(TriangleMeshFormat_plugin *plugin);
      TriangleMeshFormat_plugin *getPlugin(const char *signature) const;
      bool hasPlugin(const char *signature) const;

      static TriangleMeshFormatManager *getUniqueInstance();
    };

    /// Load and save global methods

    /// load a mesh
    LIBSL_DLL TriangleMesh             *loadTriangleMesh(const char *);

    /// load a mesh and cast it to a given format
    template<class T_VertexData,typename T_VertexFormat>
    TriangleMesh_generic<T_VertexData>* loadTriangleMesh(const char *fname)
    {
      TriangleMesh                         *mesh = loadTriangleMesh(fname);
      // check mesh has a mvf
      if (mesh->mvf().isNull()) {
        delete (mesh); // delete mesh!
        throw Errors::Fatal("loadTriangleMesh<> - No vertex format in loaded mesh! (%s)",fname);
        return NULL;
      }
      // dynamic cast
      try {
        TriangleMesh_generic<T_VertexData> *ptr  = mesh->dynamicCast< T_VertexData > ( MVF::make<T_VertexFormat>() , mesh->mvf().raw() );
        return (ptr);
      } catch (Errors::Fatal& f) {
        delete (mesh); // delete mesh!
        throw Errors::Fatal("loadTriangleMesh<> - Error during cast! (%s)",f.message());
        return NULL;
      }
    }

    /// save a mesh
    LIBSL_DLL void              saveTriangleMesh(const char *,const TriangleMesh *);

    // --------------------------------------------------------------

  } //namespace LibSL::Mesh
} //namespace LibSL

// ------------------------------------------------------

#define MESH_FORMAT_MANAGER (*LibSL::Mesh::TriangleMeshFormatManager::getUniqueInstance())

// --------------------------------------------------------------
