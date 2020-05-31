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

#include "Mesh.h"

#include <LibSL/Errors/Errors.h>
using namespace LibSL::Errors;
#include <LibSL/Memory/Array.h>
using namespace LibSL::Memory::Array;
#include <LibSL/Memory/Pointer.h>
using namespace LibSL::Memory::Pointer;
#include <LibSL/CppHelpers/CppHelpers.h>
using namespace LibSL::CppHelpers;
#include <LibSL/Geometry/AAB.h>
using namespace LibSL::Geometry;
#include <LibSL/Math/Math.h>
#include <LibSL/Math/Matrix4x4.h>
using namespace LibSL::Math;

//---------------------------------------------------------------------------

#include <algorithm>
#include <cfloat>
#include <cctype>
#include <set>
#include <map>
#include <queue>
#include <cstring>

using namespace std;

//---------------------------------------------------------------------------

#define NAMESPACE LibSL::Mesh

//---------------------------------------------------------------------------

// Mesh format manager unique instance
NAMESPACE::TriangleMeshFormatManager *NAMESPACE::TriangleMeshFormatManager::s_Manager=NULL;

//---------------------------------------------------------------------------

NAMESPACE::TriangleMeshFormatManager *NAMESPACE::TriangleMeshFormatManager::getUniqueInstance()
{
  if (s_Manager == NULL) {
    s_Manager = new TriangleMeshFormatManager();
  }
  return (s_Manager);
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMeshFormatManager::TriangleMeshFormatManager()
{

}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMeshFormatManager::registerPlugin(TriangleMeshFormat_plugin *plugin)
{
  std::map<std::string,TriangleMeshFormat_plugin *>::iterator
    P=m_Plugins.find(std::string(plugin->signature()));
  if (P == m_Plugins.end()) {
    m_Plugins[std::string(plugin->signature())]=plugin;
    // std::cerr << sprint("[TriangleMeshFormatManager] registered %s format.\n",plugin->signature());
  } else {
    // LIBSL_FATAL_ERROR_WITH_ARGS("TriangleMeshFormatManager::registerPlugin - plugin '%s' is already present !",plugin->signature());
  }
}

//---------------------------------------------------------------------------

template<class charT> charT toLower(charT c) {
  return tolower(c); // explicitely call one argument version of tolower
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMeshFormat_plugin *
NAMESPACE::TriangleMeshFormatManager::getPlugin(const char *signature) const
{
  std::string ext = std::string(signature);
  std::transform(ext.begin(),ext.end(), ext.begin(), toLower<char>);
   std::map<std::string,NAMESPACE::TriangleMeshFormat_plugin *>::const_iterator
    P=m_Plugins.find(ext);
  if (P == m_Plugins.end()) {
    LIBSL_FATAL_ERROR_WITH_ARGS("TriangleMeshFormatManager - Cannot find any plugin for '%s' (unknown format)",signature);
  } else {
    return ((*P).second);
  }
}

//---------------------------------------------------------------------------

bool NAMESPACE::TriangleMeshFormatManager::hasPlugin(const char *signature) const
{
  std::string ext = std::string(signature);
  std::transform(ext.begin(), ext.end(), ext.begin(), toLower<char>);
  std::map<std::string, NAMESPACE::TriangleMeshFormat_plugin *>::const_iterator
    P = m_Plugins.find(ext);
  if (P == m_Plugins.end()) {
    return false;
  } else {
    return true;
  }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::loadTriangleMesh(const char *fname)
{
  NAMESPACE::TriangleMeshFormatManager&
    manager=(*NAMESPACE::TriangleMeshFormatManager::getUniqueInstance());
  const char *pos=strrchr(fname,'.');
  if (pos == NULL) {
    LIBSL_FATAL_ERROR_WITH_ARGS("TriangleMesh - Cannot determine file type ('%s')",fname);
  }
  return (manager.getPlugin(pos+1)->load(fname));
}

//---------------------------------------------------------------------------

void NAMESPACE::saveTriangleMesh(const char *fname,const NAMESPACE::TriangleMesh *mesh)
{
  NAMESPACE::TriangleMeshFormatManager&
    manager=(*NAMESPACE::TriangleMeshFormatManager::getUniqueInstance());
  const char *pos=strrchr(fname,'.');
  if (pos == NULL) {
    LIBSL_FATAL_ERROR_WITH_ARGS("TriangleMesh - Cannot determine file type ('%s')",fname);
  }
  manager.getPlugin(pos+1)->save(fname,mesh);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

AABox NAMESPACE::TriangleMesh::computeBBox()
{
  float minx=1e20f,maxx=-1e20f;
  float miny=1e20f,maxy=-1e20f;
  float minz=1e20f,maxz=-1e20f;

  ForIndex(p,numVertices()) {
    float x = posAt(p)[0];
    float y = posAt(p)[1];
    float z = posAt(p)[2];
    if (x > maxx) maxx = x;
    if (x < minx) minx = x;
    if (y > maxy) maxy = y;
    if (y < miny) miny = y;
    if (z > maxz) maxz = z;
    if (z < minz) minz = z;
  }
  AAB<3> bx(v3f(minx,miny,minz), v3f(maxx,maxy,maxz));
  m_BBox         = bx;
  m_BBoxComputed = true;
  return (bx);
}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::scaleToUnitCube(float scale,bool center)
{
  AABox bx  = bbox();
  float maxd = (1.0f/scale) *
    Math::max(bx.maxCorner()[0]-bx.minCorner()[0],
    Math::max(bx.maxCorner()[1]-bx.minCorner()[1],bx.maxCorner()[2]-bx.minCorner()[2]));
  ForIndex(p,numVertices()) {
    float x = posAt(p)[0];
    float y = posAt(p)[1];
    float z = posAt(p)[2];
    posAt(p)[0] = ((x-bx.minCorner()[0])/maxd);
    posAt(p)[1] = ((y-bx.minCorner()[1])/maxd);
    posAt(p)[2] = ((z-bx.minCorner()[2])/maxd);
  }
  if (center) {
    v3 center = LibSL::Math::Triple(
      (((bx.maxCorner()[0]+bx.minCorner()[0])/2.0f-bx.minCorner()[0])/maxd),
      (((bx.maxCorner()[1]+bx.minCorner()[1])/2.0f-bx.minCorner()[1])/maxd),
      (((bx.maxCorner()[2]+bx.minCorner()[2])/2.0f-bx.minCorner()[2])/maxd));

    ForIndex(p,numVertices()) {
      posAt(p)[0] = posAt(p)[0]+0.5f-center[0];
      posAt(p)[1] = posAt(p)[1]+0.5f-center[1];
      posAt(p)[2] = posAt(p)[2]+0.5f-center[2];
    }
  }
  m_BBoxComputed = false;
}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::applyTransform(const m4x4f& trsf)
{
  ForIndex(p,numVertices()) {
    posAt(p) = trsf.mulPoint( posAt(p) );
  }
  m_BBoxComputed = false;
}

//---------------------------------------------------------------------------

v3f NAMESPACE::TriangleMesh::centerOn(const v3f& ctr)
{
  v3f   trl = ctr - bbox().center();
  ForIndex(p,numVertices()) {
    posAt(p) = posAt(p) + trl;
  }
  m_BBoxComputed = false;
  return trl;
}

//---------------------------------------------------------------------------

const AABox& NAMESPACE::TriangleMesh::bbox()
{
  if (!m_BBoxComputed) {
    computeBBox();
  }
  return (m_BBox);
}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::mergeVertices(float radius,uint gridsz)
{
  int        numv = (int)numVertices();

  // sort along longest bbox direction
  v3 ex   = bbox().extent();
  v3 dM   = V3F(1,0,0);
  float M = tupleMax(ex);
  if (ex[0] == M) dM = V3F(1,0,0);
  if (ex[1] == M) dM = V3F(0,1,0);
  if (ex[2] == M) dM = V3F(0,0,1);
  vector< pair<float,int> > allverts;
  ForIndex(v,numv) {
    allverts.push_back( make_pair( dot(posAt(v),dM) , v ) );
  }
  sort(allverts.begin(),allverts.end());

  Array<uint> merged(numv);
  ForIndex(n,numv) {
    merged[n] = n;
  }

  // merge points close to each others
  int left = 0;
  int right = 0;
  float avg_wnd_sz = 0;
  Console::progressTextInit(numv);
  ForIndex(v,numv) {
    Console::progressTextUpdate();
    int pA    = allverts[v].second;
    v3  posA  = posAt(merged[pA]);
    // advance left of comparison window
    while (left < numv-1 && dot(posAt( allverts[left].second ),dM) < dot(posA,dM)-radius) {
      left ++;
    }
    // advance right of comparison window
    while (right< numv-1 && dot(posAt( allverts[right].second ),dM) < dot(posA,dM)+radius) {
      right ++;
    }
    avg_wnd_sz += (right - left) + 1;
    // cerr << left << ',' << right << endl;
    // gather all in secondary window
    // for each within window, test
    ForRange(i,left,right) {
      int idx = allverts[i].second;
      float dist = sqLength(posA - posAt(merged[idx]));
      if (dist < radius*radius) {
        // always keep the first vertex found (hence the min)
        merged[idx] = Math::min(pA,merged[idx]);
      }
    }
  }
  Console::progressTextEnd();

  // gather used vertices
  set<uint>   used;
  ForIndex(n,numv) {
    used.insert(merged[n]);
  }

  // reorder vertices and truncate vertex array
  Array<uint> order(uint(used.size()));
  Array<uint> rank (numv);
  rank.fill(-1);
  uint n = 0;
  for(set<uint>::iterator I=used.begin() ; I!=used.end() ; I++) {
    order[n]    = (*I);
    rank [(*I)] = n;
    n ++;
  }
  reorderVerticesAndTruncate(order);

  // update triangles
  ForIndex(t,numTriangles()) {
    t_Triangle& tri = triangleAt(t);
    ForIndex(c,3) {
      tri[c]  = rank[merged[tri[c]]];
      sl_assert(tri[c] < numVertices());
    }
  }

}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::mergeVerticesExact()
{
  int        numv = (int)numVertices();

  // sort along longest bbox direction
  vector< pair<v3f, int> > allverts;
  ForIndex(v, numv) {
    allverts.push_back(make_pair(posAt(v), v));
  }
  sort(allverts.begin(), allverts.end());

  // traverse, compute new indices and store in compacted array
  Array<uint> newindex(numv);
  v3f prev = v3f(std::numeric_limits<float>::infinity());
  int i = -1;
  ForIndex(v, numv) {
    if (allverts[v].first != prev) {
      // cerr << allverts[v].first << "=?=" << prev << endl;
      i++;
    }
    newindex[allverts[v].second] = i;
    prev = allverts[v].first;
  }

  // reorder vertices and truncate vertex array
  Array<uint> order(i+1);
  ForIndex(v, numv) {
    order[newindex[v]] = v;
  }
  reorderVerticesAndTruncate(order);

  cerr << "before: " << numv << " after: " << i+1 << endl;

  // update triangles
  ForIndex(t, numTriangles()) {
    t_Triangle& tri = triangleAt(t);
    ForIndex(c, 3) {
      tri[c] = newindex[tri[c]];
      sl_assert(tri[c] < numVertices());
    }
  }

}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::swapAxes(const char swizzle[3])
{
  ForIndex(n,numVertices()) {
    v3f p = posAt(n);
    v3f s = 0;
    ForIndex(a,3) {
      int dest = 0;
      int sign = 1;
      if (swizzle[a] >= 'X' && swizzle[a] <= 'Z') {
        sign = -1;
        dest = swizzle[a] - 'X';
      } else if (swizzle[a] >= 'x' && swizzle[a] <= 'z') {
        sign = 1;
        dest = swizzle[a] - 'x';
      } else {
        sl_assert(false); // swizzle error
      }
      s[dest] = sign * p[a];
    }
    posAt(n) = s;
  }
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::TriangleMesh::clone() const
{
  TriangleMesh *newmesh = newInstance();
  newmesh->allocate(numVertices(),numTriangles());
  ForIndex(v,numVertices()) {
    memcpy(newmesh->vertexDataAt(v),vertexDataAt(v),sizeOfVertexData());
  }
  ForIndex(t,numTriangles()) {
    newmesh->triangleAt(t) = triangleAt(t);
  }
  return (newmesh);
}

//---------------------------------------------------------------------------

NAMESPACE::TriangleMesh *NAMESPACE::TriangleMesh::splitTriangles() const
{
  TriangleMesh *newmesh = newInstance();
  newmesh->allocate(numTriangles()*3,numTriangles());
  ForIndex(t,numTriangles()) {
    ForIndex(p,3) {
      newmesh->triangleAt(t)[p] = t*3+p;
      memcpy(newmesh->vertexDataAt(t*3+p),vertexDataAt(triangleAt(t)[p]),sizeOfVertexData());
    }
  }
  return (newmesh);
}

//---------------------------------------------------------------------------

void NAMESPACE::TriangleMesh::reorientTriangles()
{
  if (numTriangles() == 0) return;
  Array<bool> visited(numTriangles());
  visited.fill(false);

  // build edge info
  map<v2i, vector<int> > edge_to_tris;
  ForIndex(t, numTriangles()) {
    v3u tri = triangleAt(t);
    ForIndex(i, 3) {
      v2i e = v2i(min(tri[i], tri[(i + 1) % 3]), max(tri[i], tri[(i + 1) % 3]));
      edge_to_tris[e].push_back(t);
    }
  }

  while (true) {

    // find a lowest triangle not visited
    int next = -1;
    float minz = FLT_MAX;
    ForIndex(t,numTriangles()) {
      if (!visited[t]) {
        ForIndex(i, 3) {
          float z = posAt(triangleAt(t)[i])[2];
          if (z < minz) {
            minz = z;
            next = t;
          }
        }
      }
    }
    if (next == -1) {
      break; // done!
    }

    // orient 'next'
    v3u tn = triangleAt(next);
    v3f nrm = cross(posAt(tn[1]) - posAt(tn[0]), posAt(tn[2]) - posAt(tn[0]));
    if (nrm[2] > 0) {
      swap(tn[0], tn[2]);
      triangleAt(next) = tn;
    }

    // orient by local growth
    std::queue<int> q;
    q.push(next);
    visited[next] = true;

    while (!q.empty()) {
      // pop current
      int curid = q.front();
      v3u cur = triangleAt(curid);
      q.pop();
      // get neighbors
      ForIndex(i, 3) {
        v2i e = v2i(min(cur[i], cur[(i + 1) % 3]), max(cur[i], cur[(i + 1) % 3]));
        const vector<int>& neighs = edge_to_tris[e];
        if (neighs.size() == 2) { // only trust two-manifold edges
          ForIndex(n, neighs.size()) {
            if (!visited[neighs[n]]) {
              visited[neighs[n]] = true;
              // check orientation
              v3u tri_n = triangleAt(neighs[n]);
              bool fliped = false;
              v2i e_cur = v2i(cur[i], cur[(i + 1) % 3]);
              ForIndex(j, 3) {
                v2i e_j = v2i(tri_n[j], tri_n[(j + 1) % 3]);
                if (e_j == e_cur) {
                  fliped = true; break;
                }
              }
              if (fliped) {
                // flip
                swap(tri_n[0], tri_n[2]);
                triangleAt(neighs[n]) = tri_n;
              }
              q.push(neighs[n]);
            }
          }
        }
      }
    }
  }
}

//---------------------------------------------------------------------------
