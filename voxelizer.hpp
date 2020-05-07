
#include <algorithm>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace Voxelizer {
enum VoxelFill { Inside, Robust, None };

struct SimpleMesh {
  xt::xarray<float> points{};
  xt::xarray<int> indices{};
  xt::xtensor_fixed<float, xt::xshape<2, 3>> bounds{};

  int numVertices() const { return points.shape()[0]; }
  int numTriangles() const { return indices.shape()[0]; };

  xt::xarray<float> extent() const {
    return xt::view(bounds, 1, xt::all()) - xt::view(bounds, 0, xt::all());
  }
};

xt::xarray<float> voxelize_mesh(const SimpleMesh &mesh,
                                xt::xarray<int> voxel_resolution = {256, 256,
                                                                    256},
                                VoxelFill voxel_fill = VoxelFill::None);

xt::xarray<float> voxelize_stl(const std::string filename,
                               xt::xarray<int> voxel_resolution = {256, 256,
                                                                   256},
                               xt::xarray<float> bounds = {},
                               VoxelFill voxel_fill = VoxelFill::None);
} // namespace Voxelizer