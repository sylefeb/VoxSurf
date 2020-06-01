import pyvoxsurf
from mayavi import mlab

volume1 = pyvoxsurf.voxelize_stl("model.stl",200,[],"Robust")
print(volume1.shape)

# Visualize voxelized model
from tvtk.util.ctf import PiecewiseFunction
mlab.figure(size=(800,800))
vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(volume1))
mlab.title('Voxelized model',height=0.9,size=0.5)
mlab.orientation_axes()
otf = PiecewiseFunction()
otf.add_point(0,0)
otf.add_point(0.001, 1)
otf.add_point(1,1)
vol._otf = otf
vol._volume_property.set_scalar_opacity(otf)
mlab.show()


import pyvoxsurf
import trimesh
import numpy as np
from mayavi import mlab

mesh = trimesh.load("model.stl") # Load stl file

# Find the max and min coordinates of the mesh to form a bounding box
mesh_min_corner = [np.min(mesh.vertices[:,0]), np.min(mesh.vertices[:,1]), np.min(mesh.vertices[:,2])]
mesh_max_corner = [np.max(mesh.vertices[:,0]), np.max(mesh.vertices[:,1]), np.max(mesh.vertices[:,2])]
bounds = np.stack((mesh_min_corner,mesh_max_corner))

volume2 = pyvoxsurf.voxelize(mesh.vertices,mesh.faces,bounds,100,"Inside")
print(volume2.shape)

# Visualize voxelized model
from tvtk.util.ctf import PiecewiseFunction
mlab.figure(size=(800,800))
vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(volume2))
mlab.title('Voxelized model',height=0.9,size=0.5)
mlab.orientation_axes()
otf = PiecewiseFunction()
otf.add_point(0,0)
otf.add_point(0.001, 1)
otf.add_point(1,1)
vol._otf = otf
vol._volume_property.set_scalar_opacity(otf)
mlab.show()
