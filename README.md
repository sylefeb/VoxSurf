# VoxSurf
A simple, easily hackable C++ surface voxelizer (STL=>voxels)

Takes as input a file 'model.stl' from the source directory.
Outputs a voxel file named 'out.vox' that can be read by 'MagicaVoxel' https://ephtracy.github.io/

This is meant as an introductory, easily modifiable implementation, and includes absoltutely zero optimization (e.g. SSE and  multi-core would come to mind) nor any option such as performing a conservative voxelization. For other implementations, see links below.

The basic principle is to rasterize triangles using three 2D axis aligned grids, using integer arithmetic (fixed floating point) for robust triangle interior checks.

Very simple, CPU only, no dependencies and surprisingly efficient despite the straightforward implementation. Higher resolutions could easily be reached by not storing the voxels as a 3D array of booleans (e.g. use blocking or an octree).

Here is a relatively large model voxelized at 1024^3 in ~1.5 seconds on a Core i5-3570, 3.4GHz. (Model: [Ford engine block by Ford](https://www.thingiverse.com/thing:40257), rendered in MagicaVoxel viewer).

![voxels](vox1024.jpg)

## Links
 * CUDA voxelization (GPU, fast) https://github.com/Forceflow/cuda_voxelizer
 * Michael Schwarz and Hans-Peter Seidel paper on the topic http://research.michael-schwarz.com/publ/files/vox-siga10.pdf
 * Header only voxelization in C https://github.com/karimnaaji/voxelizer
 
## Credits
The included STL model is [3D knot by chylld](https://www.thingiverse.com/thing:5506/#files)
