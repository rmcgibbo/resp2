import numpy as np
import mdtraj as md
from cffi import FFI

ffi = FFI()
ffi.cdef("void dotsphere(int density, double* points);")
ffi.cdef("void vdw_surface(double* coordinates, char* elements, int n_elements, double scale_factor, double density, double* out, int* n_out);")
C = ffi.dlopen('dotsphere.so')

n_points = 72
points = np.empty((n_points, 3))
C.dotsphere(len(points), ffi.cast('double*', points.ctypes.data))

def plot_vertices(vertices):
    """Plot a set of vertices in 3D for debugging
    
    Example
    -------
    >>> plot_vertices(dotsphere1(100))
    """
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d import Axes3D
    fig = pp.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c='k', marker='.')
    pp.show()
    

    
  
t = md.load('native.pdb')
elements = [a.element.symbol for a in t.topology.atoms]
# convert nm to A
xyz = np.ascontiguousarray(t.xyz[0]*10, dtype=np.double)

out = np.empty((10000, 3))
out.fill(-1)
i = np.array([1], dtype=np.int32)
C.vdw_surface(ffi.cast('double*', xyz.ctypes.data), ' '.join(elements), len(elements), 2.4, 1.0,
              ffi.cast('double*', out.ctypes.data),  ffi.cast('int*', i.ctypes.data))
i = i[0]
#print xyz
print out[:i]
plot_vertices(out[:i])
