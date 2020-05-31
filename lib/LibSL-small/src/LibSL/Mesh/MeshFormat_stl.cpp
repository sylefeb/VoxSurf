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
//---------------------------------------------------------------------------
#include "LibSL.precompiled.h"
//---------------------------------------------------------------------------

#include "MeshFormat_stl.h"
using namespace LibSL::Mesh;

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;
#include <LibSL/Math/Vertex.h>
using namespace LibSL::Math;
#include <LibSL/CppHelpers/BasicParser.h>
using namespace LibSL::CppHelpers;

#include <string>

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Mesh

//---------------------------------------------------------------------------

/// Declaring a global will automatically register the plugin
namespace {
  NAMESPACE::MeshFormat_stl s_Stl;  /// FIXME: this mechanism does not work with VC++
}                                   ///        see also MeshFormatManager constructor

//---------------------------------------------------------------------------

NAMESPACE::MeshFormat_stl::MeshFormat_stl()
{
  try {
    // register plugin
    MESH_FORMAT_MANAGER.registerPlugin(this);
  } catch (LibSL::Errors::Fatal& e) {
    std::cerr << e.message() << std::endl;
  }
}

//---------------------------------------------------------------------------

using namespace std;

//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::MeshFormat_stl::loadASCII(const char *fname) const
{
  LIBSL_BEGIN;
  vector<MeshFormat_stl::t_VertexData> verts;
  vector<v3u>                          tris;

  BasicParser::FileStream stream(fname);
  BasicParser::Parser<BasicParser::FileStream> parser(stream);
  const char *s = parser.readString();
  if (strcmp(s,"solid")) {
	 throw Fatal("[MeshFormat_stl::loadASCII] - file '%s' is not an ascii STL",fname);
  }
  s = parser.readString("\n");
  cerr << sprint("[stl] loading mesh '%s', file '%s'",s,fname) << endl;
  while (!parser.eof()) {
	  if (parser.eof()) break;
	  s = parser.readString();
	  if (!strcmp(s,"endsolid")) {
		  break;
	  } else if (strcmp(s,"facet")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('facet' expected, got '%s')",s);
	  }
	  s = parser.readString();
	  if (strcmp(s,"normal")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('normal' expected, got '%s')",s);
	  }
	  MeshFormat_stl::t_VertexData vd;
	  vd.nrm[0] = parser.readFloat();
	  vd.nrm[1] = parser.readFloat();
	  vd.nrm[2] = parser.readFloat();
	  s = parser.readString();
	  if (strcmp(s,"outer")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('outer' expected, got '%s')",s);
	  }
	  s = parser.readString();
	  if (strcmp(s,"loop")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('loop' expected, got '%s')",s);
	  }
	  ForIndex(v,3) {
		s = parser.readString();
		if (strcmp(s,"vertex")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('vertex' expected, got '%s')",s);
		}
		vd.pos[0] = parser.readFloat();
		vd.pos[1] = parser.readFloat();
		vd.pos[2] = parser.readFloat();
		verts.push_back( vd );
	  }
	  int id0 = (int)verts.size()-3;
	  tris.push_back(V3U(id0,id0+1,id0+2));
	  s = parser.readString();
	  if (strcmp(s,"endloop")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('endloop' expected, got '%s')",s);
	  }
	  s = parser.readString();
	  if (strcmp(s,"endfacet")) {
		  throw Fatal("[MeshFormat_stl::loadASCII] - invalid file format ('endfacet' expected, got '%s')",s);
	  }
  }
  // allocate mesh
  TriangleMesh_generic<MeshFormat_stl::t_VertexData> *mesh
    = new TriangleMesh_generic<MeshFormat_stl::t_VertexData>((uint)verts.size(), (uint)tris.size(), 0, AutoPtr<MVF>(MVF::make<MeshFormat_stl::t_VertexFormat>()));
  ForIndex(v,verts.size()) {
	  mesh->vertexAt(v)   = verts[v];
  }
  ForIndex(t,tris.size()) {
	  mesh->triangleAt(t) = tris[t];
  }
  return mesh;
  LIBSL_END;
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::MeshFormat_stl::loadBinary(const char *fname) const
{
  LIBSL_BEGIN;
  FILE *f = NULL;
  // open file
  fopen_s(&f,fname,"rb");
  if (f == NULL) {
    throw Fatal("[MeshFormat_stl::load] - cannot open file '%s'",fname);
  }
  // read first characters
  char header[80];
  fread(header,sizeof(char),5,f);
  /*if ( ! memcmp(header,"solid",5) ) {
    throw Fatal("[MeshFormat_stl::loadBinary] - file '%s' is an ascii STL - call loadASCII instead",fname);
  }*/
  // skip header (!)
  fread(header+5,sizeof(char),75,f);
  // read number of triangles
  uint nTris = 0;
  fread(&nTris,sizeof(uint),1,f);
  cerr << "[stl] number of triangles: " << nTris << endl;
  // allocate mesh
  TriangleMesh_generic<MeshFormat_stl::t_VertexData> *mesh
    = new TriangleMesh_generic<MeshFormat_stl::t_VertexData>(nTris * 3, nTris, 0, AutoPtr<MVF>(MVF::make<MeshFormat_stl::t_VertexFormat>()));
  ForIndex(t,nTris) {
    v3f nrm,pt[3];
    fread(&nrm[0],sizeof(float),3,f);
    mesh->triangleAt(t) = V3U(t*3+0,t*3+1,t*3+2);
    fread(&pt[0][0],sizeof(float),3*3,f);
    ForIndex(v,3) {
      mesh->vertexAt(t*3+v).pos = pt[v];
      mesh->vertexAt(t*3+v).nrm = nrm;
    }
    ushort attr = 0;
    fread(&attr,sizeof(ushort),1,f);
  }
  // done with file
  fclose(f);
  cerr << "[stl done]\n";
  return mesh;
  LIBSL_END;
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::MeshFormat_stl::load(const char *fname) const
{
	LIBSL_BEGIN;
	FILE *f = NULL;
	// open file
	fopen_s(&f,fname,"rb");
	if (f == NULL) {
		throw Fatal("[MeshFormat_stl::load] - cannot open file '%s'",fname);
	}
	// read first characters
    char header[80];
	// done with file for now

	// call proper function
  NAMESPACE::TriangleMesh *mesh = NULL;
  //skip header
  fread(header,sizeof(char),80,f);
  // read number of triangles
  uint nTris = 0;
  fread(&nTris,sizeof(uint),1,f);
  long size = nTris*50+84;//size in octet of the file: headerSize(80)+4(Number of triangle)+NumberOfTriangle*50
  fseek(f, 0L, SEEK_END);

  long sz = ftell(f);// if size == sz it's a binary file
  fclose(f);//done with file
  //cerr << endl <<"real size " << sz << endl;
  //cerr << "comp size " << size << endl;
  if(size !=  sz){
  //	if ( ! memcmp(header,"solid",5) ) { ///// TODO FIXME this test is not reliable as some binary files start with solid
       mesh = loadASCII(fname);
    } else {
       mesh = loadBinary(fname);
    }
  // compute normals
  ForIndex(t,mesh->numTriangles()) {
    v3f pts[3];
    ForIndex(i,3) {
      pts[i] = mesh->posAt( mesh->triangleAt(t)[i] );
    }
    v3f nrm = normalize_safe(cross( pts[1]-pts[0] , pts[2]-pts[0] ));
    ForIndex(i,3) {
      t_VertexData *d = (t_VertexData *)mesh->vertexDataAt( mesh->triangleAt(t)[i] );
      d->nrm = nrm;
    }
  }
  return mesh;
	LIBSL_END;
}

//---------------------------------------------------------------------------

void NAMESPACE::MeshFormat_stl::save(const char *fname,const NAMESPACE::TriangleMesh *mesh) const
{
  // TODO/FIXME This assumes a vertex format compatible with MeshFormat_stl at a raw level
  //            *very* dangerous. Use a vertex attribute accessor!

  FILE *f = NULL;
  fopen_s(&f,fname,"wb");
  if (f == NULL) {
    throw Fatal("[MeshFormat_stl::save] - cannot open file '%s' for writing",fname);
  }
  // header
  char header[80];
  memset(header,0x00,80);
  sprintf(header,"LibSL_STL_1.0");
  fwrite(header,sizeof(char),80,f);
  // num triangles
  uint nTris = mesh->numTriangles();
  fwrite(&nTris,sizeof(uint),1,f);
  // vertices
  ForIndex(t,mesh->numTriangles()) {
    v3u tri = mesh->triangleAt(t);
    v3f nrm = 0;
    v3f pt[3];
    ForIndex(v,3) {
      MeshFormat_stl::t_VertexData *d = (MeshFormat_stl::t_VertexData*)mesh->vertexDataAt(tri[v]);
      pt[v] = d->pos;
      nrm  += d->nrm;
    }
    nrm = nrm / 3.0f;
    fwrite(&nrm[0],sizeof(float),3,f);
    fwrite(&pt[0][0],sizeof(float),3*3,f);
    ushort attr = 0;
    fwrite(&attr,sizeof(ushort),1,f);
  }

  fclose(f);
}

//---------------------------------------------------------------------------
