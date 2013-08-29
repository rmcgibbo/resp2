import numpy as np
from cffi import FFI

ffi = FFI()
ffi.cdef("void c_dotsphere(int density, double* points);")
C = ffi.dlopen('dotsphere.so')

n_points = 72
points = np.empty((n_points, 3))
C.c_dotsphere(len(points), ffi.cast('double*', points.ctypes.data))

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
    
print points
plot_vertices(points)
