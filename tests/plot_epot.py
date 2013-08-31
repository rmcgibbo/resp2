from StringIO import StringIO
import numpy as np
import sys
def plot_vertices(vertices, color):
    """Plot a set of vertices in 3D for debugging

    Example
    -------
    >>> plot_vertices(dotsphere1(100))
    """
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d import Axes3D
    fig = pp.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(vertices[:, 0], vertices[:, 1], vertices[:, 2], c=color, marker='o')
    #pp.colorbar()
    pp.show()


buff = []

with open(sys.argv[1]) as f:
    while not f.readline().strip() == 'Electrostatic potentials at van der Waals shells:':
        continue

    for i in range(3):
        f.readline()

    for line in f:
        if '-------' not in line:
            buff.append(line)
        else:
            break

    data = np.loadtxt(StringIO(''.join(buff)))

plot_vertices(data[:,:3], data[:,3])

