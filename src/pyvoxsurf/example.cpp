// --------------------------------------------------------------
// SL 2018-01-02
// A simple, easily hackable CPU surface voxelizer
// MIT-license
// (c) Sylvain Lefebvre, https://github.com/sylefeb
// --------------------------------------------------------------

/*

The basic principle is to rasterize triangles using three 2D axis
aligned grids, using integer arithmetic (fixed floating point)
for robust triangle interior checks.

Very simple and quite efficient despite a straightforward implementation.
Higher resolutions could easily be reached by not storing the
voxels as a 3D array of booleans (e.g. use blocking or an octree).

For the inside fill to work properly, the mesh has to be perfectly
watertight, with exactly matching vertices between neighboring
verticies.

*/

#include "voxelizer.hpp"

#include <H5Cpp.h>
#include <iostream>

static void save_xt_volume(const std::string &output_file,
                           const xt::xarray<float> &volume,
                           const std::string &dataset_name = "voxelVolume") {
  H5::H5File file(output_file, H5F_ACC_TRUNC);
  auto shape = volume.shape();
  std::vector<hsize_t> fdim = {shape[0], shape[1], shape[2]};
  std::vector<hsize_t> start = {0, 0, 0};
  std::vector<hsize_t> count = {shape[0], shape[1], shape[2]};
  std::vector<hsize_t> chunks = {std::min(32ul, volume.shape()[0]),
                                 std::min(32ul, volume.shape()[1]),
                                 std::min(32ul, volume.shape()[2])};

  float fillvalue = NAN;
  H5::DSetCreatPropList proplist;
  proplist.setDeflate(4);
  proplist.setFillValue(H5::PredType::NATIVE_FLOAT, &fillvalue);
  proplist.setChunk(chunks.size(), chunks.data());

  H5::DataSpace fspace(fdim.size(), fdim.data());
  H5::DataSet dataset = file.createDataSet(
      dataset_name, H5::PredType::NATIVE_FLOAT, fspace, proplist);
  dataset.write(volume.data(), H5::PredType::NATIVE_FLOAT);
}

int main(int argc, char **argv) {
  xt::xarray<float> voxels;
  if (argc >= 2) {
    voxels = Voxelizer::voxelize_stl(argv[1], 128);
  } else {
    Voxelizer::SimpleMesh mesh;
    xt::xarray<float> cube_verts = {{-1.1, 1.1, 1.1},  {1.1, -1.1, 1.1},
                                    {1.1, 1.1, 1.1},   {-1.1, -1.1, -1.1},
                                    {1.1, -1.1, -1.1}, {-1.1, -1.1, 1.1},
                                    {-1.1, 1.1, -1.1}, {1.1, 1.1, -1.1}};
    mesh.points = cube_verts;
    mesh.indices = {{0, 1, 2}, {1, 3, 4}, {5, 6, 3}, {7, 3, 6},
                    {2, 4, 7}, {0, 7, 6}, {0, 5, 1}, {1, 5, 3},
                    {5, 0, 6}, {7, 4, 3}, {2, 1, 4}, {0, 2, 7}};

    mesh.bounds = {{-1.51, -1.51, -1.51}, {1.51, 1.51, 1.51}};
    voxels = Voxelizer::voxelize_mesh(mesh, 512);
  }
  save_xt_volume("out.h5", voxels);
}
