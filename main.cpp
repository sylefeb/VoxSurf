// --------------------------------------------------------------
// SL 2018-01-02
// A simple, easily hackable CPU surface voxelizer
// MIT-license
// (c) Sylvain Lefebvre, https://github.com/sylefeb
/*

Takes as input a file 'model.stl' from the source directory.
Outputs a voxel file named 'out.vox' that can be read by 'MagicaVoxel' https://ephtracy.github.io/

Change VOXEL_RESOLUTION to fit your needs.
Set    VOXEL_INSIDE to one to fill in the interior

The basic principle is to rasterize triangles using three 2D axis
aligned grids, using integer arithmetic (fixed floating point)
for robust triangle interior checks.

Very simple and quite efficient despite a straightforward implementation.
Higher resolutions could easily be reached by not storing the
voxels as a 3D array of booleans (e.g. use blocking or an octree).

*/

#include <LibSL/LibSL.h>

#include <iostream>
#include <algorithm>
#include <queue>
using namespace std;

#include "path.h"

// --------------------------------------------------------------

#define VOXEL_RESOLUTION 128
#define VOXEL_INSIDE     0

// --------------------------------------------------------------

#define FP_POW    16
#define FP_SCALE  (1<<FP_POW)
#define BOX_SCALE v3f(VOXEL_RESOLUTION*FP_SCALE)

// --------------------------------------------------------------

// saves a voxel file (.vox format, can be imported by MagicaVoxel)
void saveAsVox(const char *fname, const Array3D<bool>& voxs)
{
  Array<v3b> palette(256); // RGB palette
  palette.fill(0);
  palette[127] = v3b(255, 255, 255);
  FILE *f;
  f = fopen(fname, "wb");
  sl_assert(f != NULL);
  long sx = voxs.xsize(), sy = voxs.ysize(), sz = voxs.zsize();
  fwrite(&sx, 4, 1, f);
  fwrite(&sy, 4, 1, f);
  fwrite(&sz, 4, 1, f);
  ForIndex(i, sx) {
    ForIndex(j, sy) {
      ForRangeReverse(k, sz - 1, 0) {
        uchar pal = voxs.at(i, j, k) ? 127 : 255;
        fwrite(&pal, sizeof(uchar), 1, f);
      }
    }
  }
  fwrite(palette.raw(), sizeof(v3b), 256, f);
  fclose(f);
}

// --------------------------------------------------------------

inline bool isInTriangle(int i, int j, const v3i& p0, const v3i& p1, const v3i& p2, int& _depth)
{
  v2i delta_p0 = v2i(i, j) - v2i(p0);
  v2i delta_p1 = v2i(i, j) - v2i(p1);
  v2i delta_p2 = v2i(i, j) - v2i(p2);
  v2i delta10 = v2i(p1) - v2i(p0);
  v2i delta21 = v2i(p2) - v2i(p1);
  v2i delta02 = v2i(p0) - v2i(p2);

  int64_t c0 = (int64_t)delta_p0[0] * (int64_t)delta10[1] - (int64_t)delta_p0[1] * (int64_t)delta10[0];
  int64_t c1 = (int64_t)delta_p1[0] * (int64_t)delta21[1] - (int64_t)delta_p1[1] * (int64_t)delta21[0];
  int64_t c2 = (int64_t)delta_p2[0] * (int64_t)delta02[1] - (int64_t)delta_p2[1] * (int64_t)delta02[0];
  bool inside = (c0 <= 0 && c1 <= 0 && c2 <= 0) || (c0 >= 0 && c1 >= 0 && c2 >= 0);

  if (inside) {
    int64_t area = c0 + c1 + c2;
    int64_t b0 = (c1 << 10) / area;
    int64_t b1 = (c2 << 10) / area;
    int64_t b2 = (1 << 10) - b0 - b1;
    _depth = (int)((b0 * p0[2] + b1 * p1[2] + b2 * p2[2]) >> 10);
  }
  return inside;
}

// --------------------------------------------------------------

class swizzle_xyz
{
public:
  inline v3i forward(const v3i& v)  const { return v; }
  inline v3i backward(const v3i& v) const { return v; }
};

class swizzle_zxy
{
public:
  inline v3i forward(const v3i& v)  const { return v3i(v[2], v[0], v[1]); }
  inline v3i backward(const v3i& v) const { return v3i(v[1], v[2], v[0]); }
};

class swizzle_yzx
{
public:
  inline v3i forward(const v3i& v)  const { return v3i(v[1], v[2], v[0]); }
  inline v3i backward(const v3i& v) const { return v3i(v[2], v[0], v[1]); }
};

// --------------------------------------------------------------

template <class S>
void rasterize(
  const v3u&                  tri,
  const std::vector<v3i>&     pts,
  Array3D<bool>&             _voxs)
{
  const S swizzler;
  v3i tripts[3] = {
    swizzler.forward(pts[tri[0]]),
    swizzler.forward(pts[tri[1]]),
    swizzler.forward(pts[tri[2]])
  };
  // check if triangle is valid
  v2i delta10 = v2i(tripts[1]) - v2i(tripts[0]);
  v2i delta21 = v2i(tripts[2]) - v2i(tripts[1]);
  v2i delta02 = v2i(tripts[0]) - v2i(tripts[2]);
  if (delta10 == v2i(0)) return;
  if (delta21 == v2i(0)) return;
  if (delta02 == v2i(0)) return;
  if (delta02[0] * delta10[1] - delta02[1] * delta10[0] == 0) return;
  // proceed
  AAB<2, int> pixbx;
  pixbx.addPoint(v2i(tripts[0]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[1]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[2]) / FP_SCALE);
  for (int j = pixbx.minCorner()[1]; j <= pixbx.maxCorner()[1]; j++) {
    for (int i = pixbx.minCorner()[0]; i <= pixbx.maxCorner()[0]; i++) {
      int depth;
      if (isInTriangle(
        (i << FP_POW) + (1 << (FP_POW - 1)), // centered
        (j << FP_POW) + (1 << (FP_POW - 1)), // centered
        tripts[0], tripts[1], tripts[2], depth)) {
        v3i vx = swizzler.backward(v3i(i, j, depth >> FP_POW));
        // tag the voxel as occupied
        // NOTE: voxels are likely to be tagged multiple times (e.g. center exactly on edge, overlaps, etc.)
        _voxs.at(vx[0], vx[1], vx[2]) = true;
      }
    }
  }
}

// --------------------------------------------------------------

void fillInside(Array3D<bool>& _voxs)
{
  const v3i neighs[6] = { v3i(1,0,0),v3i(-1,0,0),v3i(0,1,0),v3i(0,-1,0),v3i(0,0,1),v3i(0,0,-1) };
  // track voxels that are outside
  Array3D<bool> outside(VOXEL_RESOLUTION, VOXEL_RESOLUTION, VOXEL_RESOLUTION);
  outside.fill(false);
  // assumes origin is outside
  std::queue<v3i> q;
  q.push(v3i(0, 0, 0));
  while (!q.empty()) {
    // pop next
    v3i cur = q.front();
    q.pop();
    // tag as outside
    outside.at(cur[0], cur[1], cur[2]) = true;
    // propagate
    ForIndex(i, 6) {
      v3i n = cur + neighs[i];
      if ( n[0] >= 0 && n[0] < VOXEL_RESOLUTION
        && n[1] >= 0 && n[1] < VOXEL_RESOLUTION 
        && n[2] >= 0 && n[2] < VOXEL_RESOLUTION) {
        if (!_voxs.at(n[0], n[1], n[2]) && !outside.at(n[0], n[1], n[2])) {
          outside.at(n[0], n[1], n[2]) = true;
          q.push(n);
        }
      }
    }
  }
  // set inner voxels
  ForArray3D(_voxs, i, j, k) {
    if (!_voxs.at(i, j, k)) {
      _voxs.at(i, j, k) = !outside.at(i, j, k);
    }
  }
}

// --------------------------------------------------------------

int main(int argc, char **argv)
{

  try {

    // load triangle mesh
    TriangleMesh_Ptr mesh(loadTriangleMesh(SRC_PATH "/model.stl"));    
    // produce (fixed fp) integer vertices and triangles
    std::vector<v3i> pts;
    std::vector<v3u> tris;
    {
      m4x4f boxtrsf = scaleMatrix(BOX_SCALE) * translationMatrix(v3f(0.5f))
        * scaleMatrix(v3f(0.95f) / tupleMax(mesh->bbox().extent()))
        * translationMatrix(-mesh->bbox().center());
      // transform vertices
      pts.resize(mesh->numVertices());
      ForIndex(p, mesh->numVertices()) {
        v3f pt   = mesh->posAt(p);
        v3f bxpt = boxtrsf.mulPoint(pt);
        v3i ipt  = v3i(clamp(round(bxpt), v3f(0.0f), BOX_SCALE - v3f(1.0f)));
        pts[p]   = ipt;
      }
      // prepare triangles
      tris.reserve(mesh->numTriangles());
      ForIndex(t, mesh->numTriangles()) {
        v3u tri = mesh->triangleAt(t);
        tris.push_back(tri);
      }
    }
    
    // rasterize into voxels
    Array3D<bool> voxs(VOXEL_RESOLUTION, VOXEL_RESOLUTION, VOXEL_RESOLUTION);
    voxs.fill(false);
    {
      Timer tm("rasterization");
      Console::progressTextInit((int)tris.size());
      ForIndex(t, tris.size()) {
        Console::progressTextUpdate(t);
        rasterize<swizzle_xyz>(tris[t], pts, voxs); // xy view
        rasterize<swizzle_yzx>(tris[t], pts, voxs); // yz view
        rasterize<swizzle_zxy>(tris[t], pts, voxs); // zx view
      }
      Console::progressTextEnd();
      cerr << endl;
    }

    // add inner voxels
#if VOXEL_INSIDE
    {
      Timer tm("flood fill");
      cerr << "filling in/out ... ";
      fillInside(voxs);
      cerr << " done." << endl;
    }
#endif

    // save the result
    saveAsVox(SRC_PATH "/out.vox", voxs);

    // report some stats
    int num_in_vox = 0;
    ForArray3D(voxs, i, j, k) {
      if (voxs.at(i, j, k)) {
        num_in_vox++;
      }
    }
    cerr << "number of set voxels: " << num_in_vox << endl;

  } catch (Fatal& e) {
    cerr << "[ERROR] " << e.message() << endl;
  }

}

/* -------------------------------------------------------- */
