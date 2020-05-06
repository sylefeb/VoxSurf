#include "voxelizer.hpp"
#include <LibSL/LibSL.h>
#include <xtensor/xindex_view.hpp>

using namespace std;

namespace Voxelizer {
static constexpr int32_t FP_POW = 16;
static constexpr int32_t FP_SCALE = 1 << FP_POW;

static constexpr uchar ALONG_X = 1;
static constexpr uchar ALONG_Y = 2;
static constexpr uchar ALONG_Z = 4;
static constexpr uchar INSIDE = 8;
static constexpr uchar INSIDE_X = 16;
static constexpr uchar INSIDE_Y = 32;
static constexpr uchar INSIDE_Z = 64;

inline bool isInTriangle(int i, int j, const v3i &p0, const v3i &p1,
                         const v3i &p2, int &_depth) {
  v2i delta_p0 = v2i(i, j) - v2i(p0);
  v2i delta_p1 = v2i(i, j) - v2i(p1);
  v2i delta_p2 = v2i(i, j) - v2i(p2);
  v2i delta10 = v2i(p1) - v2i(p0);
  v2i delta21 = v2i(p2) - v2i(p1);
  v2i delta02 = v2i(p0) - v2i(p2);

  int64_t c0 = (int64_t)delta_p0[0] * (int64_t)delta10[1] -
               (int64_t)delta_p0[1] * (int64_t)delta10[0];
  int64_t c1 = (int64_t)delta_p1[0] * (int64_t)delta21[1] -
               (int64_t)delta_p1[1] * (int64_t)delta21[0];
  int64_t c2 = (int64_t)delta_p2[0] * (int64_t)delta02[1] -
               (int64_t)delta_p2[1] * (int64_t)delta02[0];
  bool inside =
      (c0 <= 0 && c1 <= 0 && c2 <= 0) || (c0 >= 0 && c1 >= 0 && c2 >= 0);

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

class swizzle_xyz {
public:
  inline v3i forward(const v3i &v) const { return v; }
  inline v3i backward(const v3i &v) const { return v; }
  inline int along() const { return ALONG_Z; }
};

class swizzle_zxy {
public:
  inline v3i forward(const v3i &v) const { return v3i(v[2], v[0], v[1]); }
  inline v3i backward(const v3i &v) const { return v3i(v[1], v[2], v[0]); }
  inline uchar along() const { return ALONG_Y; }
};

class swizzle_yzx {
public:
  inline v3i forward(const v3i &v) const { return v3i(v[1], v[2], v[0]); }
  inline v3i backward(const v3i &v) const { return v3i(v[2], v[0], v[1]); }
  inline uchar along() const { return ALONG_X; }
};

template <class S>
static void rasterize(const v3u &tri, const std::vector<v3i> &pts,
                      Array3D<uchar> &_voxs) {
  const S swizzler;
  v3i tripts[3] = {swizzler.forward(pts[tri[0]]), swizzler.forward(pts[tri[1]]),
                   swizzler.forward(pts[tri[2]])};
  // check if triangle is valid
  v2i delta10 = v2i(tripts[1]) - v2i(tripts[0]);
  v2i delta21 = v2i(tripts[2]) - v2i(tripts[1]);
  v2i delta02 = v2i(tripts[0]) - v2i(tripts[2]);
  if (delta10 == v2i(0))
    return;
  if (delta21 == v2i(0))
    return;
  if (delta02 == v2i(0))
    return;
  if (delta02[0] * delta10[1] - delta02[1] * delta10[0] == 0)
    return;
  // proceed
  AAB<2, int> pixbx;
  pixbx.addPoint(v2i(tripts[0]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[1]) / FP_SCALE);
  pixbx.addPoint(v2i(tripts[2]) / FP_SCALE);
  for (int j = pixbx.minCorner()[1]; j <= pixbx.maxCorner()[1]; j++) {
    for (int i = pixbx.minCorner()[0]; i <= pixbx.maxCorner()[0]; i++) {
      int depth;
      if (isInTriangle((i << FP_POW) + (1 << (FP_POW - 1)), // centered
                       (j << FP_POW) + (1 << (FP_POW - 1)), // centered
                       tripts[0], tripts[1], tripts[2], depth)) {
        v3i vx = swizzler.backward(v3i(i, j, depth >> FP_POW));
        // tag the voxel as occupied
        // NOTE: voxels are likely to be hit multiple times (e.g. thin features)
        //       we flip the bit every time a hit occurs in a voxel
        _voxs.at(vx[0], vx[1], vx[2]) =
            (_voxs.at(vx[0], vx[1], vx[2]) ^ swizzler.along());
      }
    }
  }
}

// This version is more robust by using all three direction
// and voting among them to decide what is filled or not
static void fillInsideVoting(Array3D<uchar> &_voxs) {
  // along x
  for (int k = 0; k < _voxs.zsize(); k++) {
    for (int j = 0; j < _voxs.ysize(); j++) {
      bool inside = false;
      for (int i = 0; i < _voxs.xsize(); i++) {
        if (_voxs.at(i, j, k) & ALONG_X) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(i, j, k) |= INSIDE_X;
        }
      }
    }
  }
  // along y
  for (int k = 0; k < _voxs.zsize(); k++) {
    for (int j = 0; j < _voxs.xsize(); j++) {
      bool inside = false;
      for (int i = 0; i < _voxs.ysize(); i++) {
        if (_voxs.at(j, i, k) & ALONG_Y) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(j, i, k) |= INSIDE_Y;
        }
      }
    }
  }
  // along z
  for (int k = 0; k < _voxs.ysize(); k++) {
    for (int j = 0; j < _voxs.xsize(); j++) {
      bool inside = false;
      for (int i = 0; i < _voxs.zsize(); i++) {
        if (_voxs.at(j, k, i) & ALONG_Z) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(j, k, i) |= INSIDE_Z;
        }
      }
    }
  }
  // voting
  ForArray3D(_voxs, i, j, k) {
    uchar v = _voxs.at(i, j, k);
    int votes = ((v & INSIDE_X) ? 1 : 0) + ((v & INSIDE_Y) ? 1 : 0) +
                ((v & INSIDE_Z) ? 1 : 0);
    // clean
    _voxs.at(i, j, k) &= ~(INSIDE_X | INSIDE_Y | INSIDE_Z);
    if (votes > 1) {
      // tag as inside
      _voxs.at(i, j, k) |= INSIDE;
    }
  }
}

static void fillInside(Array3D<uchar> &_voxs) {
  for (int k = 0; k < _voxs.zsize(); k++) {
    for (int j = 0; j < _voxs.ysize(); j++) {
      bool inside = false;
      for (int i = 0; i < _voxs.xsize(); i++) {
        if (_voxs.at(i, j, k) & ALONG_X) {
          inside = !inside;
        }
        if (inside) {
          _voxs.at(i, j, k) |= INSIDE;
        }
      }
    }
  }
}

template <typename T> xt::xarray<T> sl2xt(Array3D<T> voxels) {
  static xt::xarray<T> xtvoxels =
      xt::zeros<T>({voxels.xsize(), voxels.ysize(), voxels.zsize()});
  auto total_size = voxels.xsize() * voxels.ysize() * voxels.zsize();
  for (size_t i = 0; i < voxels.xsize(); i++) {
    for (size_t j = 0; j < voxels.ysize(); j++) {
      for (size_t k = 0; k < voxels.zsize(); k++) {
        xtvoxels(i, j, k) = voxels.at(i, j, k);
      }
    }
  }
  return std::move(xtvoxels);
}

static xt::xarray<float> voxelize(const std::vector<v3i> &pts,
                                  const std::vector<v3u> &tris, v3u resolution,
                                  VoxelFill voxel_fill) {
  std::cout << "VOXELIZEE " << resolution << " " << pts.size() << " "
            << tris.size() << "\n ";
  Array3D<uchar> voxs(resolution);
  voxs.fill(0);
  {
    for (int t = 0; t < tris.size(); t++) {
      rasterize<swizzle_xyz>(tris[t], pts, voxs); // xy view
      rasterize<swizzle_yzx>(tris[t], pts, voxs); // yz view
      rasterize<swizzle_zxy>(tris[t], pts, voxs); // zx view
    }
  }

  // add inner voxels
  switch (voxel_fill) {
  case Inside:
    fillInside(voxs);
    break;
  case Robust:
    fillInsideVoting(voxs);
  case None:
  default:
    break;
  }
  auto xtvoxels = sl2xt(voxs);
  xt::filter(xtvoxels, xtvoxels > 0) = 1.0f;
  return std::move(xtvoxels);
}

template <typename T> auto xt2v3(xt::xarray<T> v) {
  return LibSL::Math::Tuple<T, 3>(v[0], v[1], v[2]);
}

xt::xarray<float> voxelize_mesh(const SimpleMesh &mesh,
                                uint32_t voxel_resolution,
                                VoxelFill voxel_fill) {

  // produce (fixed fp) integer vertices and triangles
  std::vector<v3i> pts;
  std::vector<v3u> tris;
  {
    const v3f box_scale = v3f(voxel_resolution * FP_SCALE);

    float factor = 0.95f;

    m4x4f boxtrsf =
        scaleMatrix(box_scale) *
        scaleMatrix(v3f(1.f) / xt::amax(mesh.extent())[0]) *
        translationMatrix((1 - factor) * 0.5f * xt2v3<float>(mesh.extent())) *
        scaleMatrix(v3f(factor)) *
        translationMatrix(xt2v3<float>(xt::view(-mesh.bounds, 0, xt::all())));

    // transform vertices
    pts.resize(mesh.numVertices());
    for (int p = 0; p < mesh.numVertices(); p++) {
      v3f pt = xt2v3<float>(xt::view(mesh.points, p, xt::all()));
      v3f bxpt = boxtrsf.mulPoint(pt);
      v3i ipt = v3i(clamp(round(bxpt), v3f(0.0f), box_scale - v3f(1.0f)));
      pts[p] = ipt;
    }
    // prepare triangles
    tris.reserve(mesh.numTriangles());
    for (int t = 0; t < mesh.numTriangles(); t++) {
      v3u tri = xt2v3<uint>(xt::view(mesh.indices, t, xt::all()));
      tris.push_back(tri);
    }
  }

  // rasterize into voxels
  v3u resolution(xt2v3<uint>(mesh.extent() / xt::amax(mesh.extent())[0] *
                             float(voxel_resolution)));

  return std::move(voxelize(pts, tris, resolution, voxel_fill));
}

xt::xarray<float> voxelize_stl(const std::string stl_file,
                               uint32_t voxel_resolution,
                               VoxelFill voxel_fill) {
  // LibSL needs to register the mesh format, but it doesn't do it
  // automatically from the library, so we create this dummy instance
  LibSL::Mesh::MeshFormat_stl _plugin_register_workaround;
  // load triangle mesh
  TriangleMesh_Ptr mesh(loadTriangleMesh(stl_file.c_str()));
  // produce (fixed fp) integer vertices and triangles
  std::vector<v3i> pts;
  std::vector<v3u> tris;
  {
    const v3f box_scale = v3f(voxel_resolution * FP_SCALE);

    float factor = 0.95f;
    m4x4f boxtrsf =
        scaleMatrix(box_scale) *
        scaleMatrix(v3f(1.f) / tupleMax(mesh->bbox().extent())) *
        translationMatrix((1 - factor) * 0.5f * mesh->bbox().extent()) *
        scaleMatrix(v3f(factor)) * translationMatrix(-mesh->bbox().minCorner());

    // transform vertices
    pts.resize(mesh->numVertices());
    for (int p = 0; p < mesh->numVertices(); p++) {
      v3f pt = mesh->posAt(p);
      v3f bxpt = boxtrsf.mulPoint(pt);
      v3i ipt = v3i(clamp(round(bxpt), v3f(0.0f), box_scale - v3f(1.0f)));
      pts[p] = ipt;
    }
    // prepare triangles
    tris.reserve(mesh->numTriangles());
    for (int t = 0; t < mesh->numTriangles(); t++) {
      v3u tri = mesh->triangleAt(t);
      tris.push_back(tri);
    }
  }

  // rasterize into voxels
  v3u resolution(mesh->bbox().extent() / tupleMax(mesh->bbox().extent()) *
                 float(voxel_resolution));

  return std::move(voxelize(pts, tris, resolution, voxel_fill));
}
} // namespace Voxelizer