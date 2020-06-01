#include "pybind11/pybind11.h"

#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pyvectorize.hpp>

#include "voxelizer.hpp"
#include <iostream>
namespace py = pybind11;
using namespace pybind11::literals;

inline Voxelizer::VoxelFill str2vfill(const std::string &voxel_fill) {
  using Voxelizer::VoxelFill;
  if (voxel_fill == "None") {
    return VoxelFill::None;
  } else if (voxel_fill == "Inside") {
    return VoxelFill::Inside;
  } else if (voxel_fill == "Robust") {
    return VoxelFill::Robust;
  }
  return VoxelFill::None;
}

inline xt::pyarray<float> voxelize_stl(const std::string &filename,
                                       const int resolution,
                                       const xt::pyarray<float> &bounds = {},
                                       const std::string &voxel_fill = "None") {
  return Voxelizer::voxelize_stl(filename, resolution, bounds,
                                 str2vfill(voxel_fill));
}

inline xt::pyarray<float> voxelize(const xt::pyarray<float> &vertices,
                                   const xt::pyarray<int> &triangle_indices,
                                   const xt::pyarray<float> &bounds,
                                   const int resolution,
                                   const std::string &voxel_fill = "None") {
  Voxelizer::SimpleMesh mesh;
  mesh.points = vertices;
  mesh.indices = triangle_indices;
  mesh.bounds = bounds;
  return Voxelizer::voxelize_mesh(mesh, resolution, str2vfill(voxel_fill));
}

// Python Module and Docstrings
PYBIND11_MODULE(pyvoxsurf, m) {
  xt::import_numpy();

  m.doc() = R"pbdoc(
        Python binding to voxelize raw meshes and stl files
    )pbdoc";

  m.def("voxelize_stl", voxelize_stl,
        "Loads an stl file and voxelizes. It will voxelize the whole file if "
        "bounds is left empty, or the specified bounds otherwise.",
        "filename"_a, "voxel_resolution"_a, "bounds"_a = xt::xarray<float>(),
        "voxel_fill"_a = "None");

  m.def("voxelize", voxelize, "Voxelizes the mesh within the specified bounds.",
        "vertices"_a, "triangle_indices"_a, "bounds"_a, "voxel_resolution"_a,
        "voxel_fill"_a = "None");
}
