# This code is adapted from GROMACS
#
# Copyright (c) 1991-2000, University of Groningen, The Netherlands.
# Copyright (c) 2001-2007, The GROMACS development team,
# check out http://www.gromacs.org for more information.
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
"""
Two routines for computing a set of dots (approximately) unfiromly
distributed on the unit sphere by repeated truncation of icosahedrons.

The two methods, dotsphere1 and dotsphere2, use slightly different procedures
which causes them to be capable of yielding different numbers of points.
"""

#############################################################################
# Imports
#############################################################################

from __future__ import division
import numpy as np

#############################################################################
# Utilities
#############################################################################


def to_rad(deg):
    "Convert degrees to radians"
    return deg * np.pi / 180.0

#############################################################################
# Globals
#############################################################################

R_H = np.sqrt(1.0 - 2.0*np.cos(to_rad(72.)))/(1.-np.cos(to_rad(72.)))
DP_TOL = 0.001

__all__ = ['dotsphere1', 'dotsphere2', 'dotsphere']

#############################################################################
# Functions
#############################################################################

def _divarc(xyz1, xyz2, div1, div2):
    """Divide an arc based on the great circle"""

    xd = xyz1[1]*xyz2[2] - xyz2[1]*xyz1[2]
    yd = xyz1[2]*xyz2[0] - xyz2[2]*xyz1[0]
    zd = xyz1[0]*xyz2[1] - xyz2[0]*xyz1[1]
    dd = np.sqrt(xd*xd + yd*yd + zd*zd)
    if dd < DP_TOL:
        raise ValueError("_divarc: rotation axis of length %f" % dd)

    d1 = sum(xyz1**2)
    if d1 < 0.5:
        raise ValueError("_divarc: vector 1 of sq.length %f" % d1)

    d2 = sum(xyz2**2)
    if d2 < 0.5:
        raise ValueError("_divarc: vector 2 of sq.length %f", d2)

    phi = np.sin(dd / np.sqrt(d1*d2))
    phi = phi * div1 / div2  # float division
    sphi = np.sin(phi)
    cphi = np.cos(phi)
    s = (xyz1[0]*xd + xyz1[0]*yd + xyz1[2]*zd) / dd

    x = xd*s*(1.0-cphi)/dd + xyz1[0] * cphi + (yd*xyz1[2] - xyz1[1]*zd)*sphi/dd
    y = yd*s*(1.0-cphi)/dd + xyz1[1] * cphi + (zd*xyz1[0] - xyz1[2]*xd)*sphi/dd
    z = zd*s*(1.0-cphi)/dd + xyz1[2] * cphi + (xd*xyz1[1] - xyz1[0]*yd)*sphi/dd
    dd = np.sqrt(x*x + y*y + z*z)

    return np.array([x/dd, y/dd, z/dd])


def icosahedron_vertices():
    """Compute the vertices of an icosahedron

    This code was adapted from GROMACS's nsc.c, distributed under the GNU
    LGPL. See this file's header for the copyright information.

    Returns
    -------
    verts : np.ndarray, shape=(12, 3)
        Cartesian coordinates of the 12 verticles of a unit icosahedron.
    """

    rg = np.cos(to_rad(72.))/(1.-np.cos(to_rad(72.)))

    verts = np.empty((12, 3))
    verts[0] = [0.0, 0.0, 1.0]
    verts[1] = [R_H*np.cos(to_rad(72.)), R_H*np.sin(to_rad(72.)), rg]
    verts[2] = [R_H*np.cos(to_rad(144.)), R_H*np.sin(to_rad(144.)), rg]
    verts[3] = [R_H*np.cos(to_rad(216.)), R_H*np.sin(to_rad(216.)), rg]
    verts[4] = [R_H*np.cos(to_rad(288.)), R_H*np.sin(to_rad(288.)), rg]
    verts[5] = [R_H, 0, rg]
    verts[6] = [R_H*np.cos(to_rad(36.)), R_H*np.sin(to_rad(36.)), -rg]
    verts[7] = [R_H*np.cos(to_rad(108.)), R_H*np.sin(to_rad(108.)), -rg]
    verts[8] = [-R_H, 0, -rg]
    verts[9] = [R_H*np.cos(to_rad(252.)), R_H*np.sin(to_rad(252.)), -rg]
    verts[10] = [R_H*np.cos(to_rad(324.)), R_H*np.sin(to_rad(324.)), -rg]
    verts[11] = [0, 0, -1]

    return verts


def dotsphere1(density):
    """Create a dot distribution over the unit shpere based on repeated
    splitting and refining the arcs of an icosahedron.

    In general, by my visual inspection (RTM, August 2013), the dot
    distributions produced by this method look worse than those produced by
    the alternative procedure, `dotsphere_icos2`. This method generally
    produces fewer points.

    Parameters
    ----------
    density : int
        Required number of dots on the unit sphere

    Returns
    -------
    dots : np.ndarray, shape=(N, 3), dtype=np.double
        Dots on the surface of the unit sphere. The number of dots will be
        at minimum equal to the `density` argument, but will be roughly two
        times larger.

    Notes
    -----
    This code was adapted from the function 'ico_dot_arc' in GROMACS's nsc.c,
    distributed under the GNU LGPL. See this file's header for the copyright
    information.

    See Also
    --------
    dotsphere_icos2 : acomplished the same goal, but based on splitting
        the faces. The two procedures are capable of yielding different
        number of points because of the different algorithms used.
    """

    # calculate tessalation level
    a = np.sqrt((density - 2.0) / 10.0)
    tess = int(np.ceil(a))

    # make vertices a LIST of numpy arrays so that we can
    # easily append to it
    vertices = [v for v in icosahedron_vertices()]

    if tess > 1:
        a = R_H*R_H*2.0*(1.0 - np.cos(to_rad(72.0)))
        # Calculate tessalation of icosahedron edges
        for i in range(11):
            for j in range(i+1, 12):
                d = sum((vertices[i] - vertices[j])**2)
                if abs(a-d) > DP_TOL:
                    continue
                for tl in range(tess):
                    vertices.append(_divarc(vertices[i], vertices[j], tl, tess))

    # Calculate tessalation of icosahedron faces
    for i in range(10):
        for j in range(i+1, 11):
            d = sum((vertices[i] - vertices[j])**2)
            if abs(a-d) > DP_TOL:
                continue

            for k in range(j+1, 12):
                d_ik = sum((vertices[i] - vertices[k])**2)
                d_jk = sum((vertices[j] - vertices[k])**2)
                if (abs(a - d_ik) > DP_TOL) or (abs(a - d_jk) > DP_TOL):
                    continue
                for tl in range(1, tess-1):
                    ji = _divarc(vertices[j], vertices[i], tl, tess)
                    ki = _divarc(vertices[k], vertices[i], tl, tess)

                    for tl2 in range(1, tess-tl):
                        ij = _divarc(vertices[i], vertices[j], tl2, tess)
                        kj = _divarc(vertices[k], vertices[j], tl2, tess)
                        ik = _divarc(vertices[i], vertices[k], tess-tl-tl2, tess)
                        jk = _divarc(vertices[j], vertices[k], tess-tl-tl2, tess)

                        xyz1 = _divarc(ki, ji, tl2, tess-tl)
                        xyz2 = _divarc(kj, ij, tl, tess-tl2)
                        xyz3 = _divarc(jk, ik, tl, tl+tl2)

                        x = xyz1 + xyz2 + xyz3
                        d = np.sqrt(sum(x**2))
                        vertices.append(x/d)

    return np.array(vertices)


def dotsphere2(density):
    """Create a dot distribution over the unit shpere based on repeated
    truncating and refining the faces of an icosahedron.

    In general, by my visual inspection (RTM, August 2013), the dot
    distributions produced by this method look "better" than those produced by
    the alternative procedure, `dotsphere_icos1`. But this method tends to
    produce more points.

    Parameters
    ----------
    density : int
        Required number of dots on the unit sphere

    Returns
    -------
    dots : np.ndarray, shape=(N, 3), dtype=np.double
        Dots on the surface of the unit sphere. The number of dots will be
        at minimum equal to the `density` argument, but will be roughly two
        times larger.

    Notes
    -----
    This code was adapted from the function 'ico_dot_dod' in GROMACS's nsc.c,
    distributed under the GNU LGPL. See this file's header for the copyright
    information.

    See Also
    --------
    dotsphere_icos1 : acomplished the same goal, but based on splitting
        the edges. The two procedures are capable of yielding different
        number of points because of the different algorithms used.
    """

    a = np.sqrt((density - 2.0) / 30.0)
    tess = max(int(np.ceil(a)), 1)

    # make vertices a LIST of numpy arrays so that we can
    # easily append to it
    vertices = [v for v in icosahedron_vertices()]

    a = R_H*R_H * 2.0 * (1.0 - np.cos(to_rad(72.0)))

    # dodecaeder vertices
    for i in range(10):
        for j in range(i+1, 11):
            d = sum((vertices[i] - vertices[j])**2)
            if abs(a-d) > DP_TOL:
                continue
            for k in range(j+1, 12):
                d_ik = sum((vertices[i] - vertices[k])**2)
                d_jk = sum((vertices[j] - vertices[k])**2)
                if (abs(a - d_ik) > DP_TOL) or (abs(a - d_jk) > DP_TOL):
                    continue

                xyz = vertices[i] + vertices[j] + vertices[k]
                d = np.sqrt(sum(xyz**2))
                vertices.append(xyz/d)

    if tess > 1:
        # square of the edge of an dodecaeder
        adod = 4.0 * (np.cos(to_rad(108.)) - np.cos(to_rad(120.))) / (1.0 - np.cos(to_rad(120.)))
        # square of the distance of two adjacent vertices of ico- and dodecaeder */
        ai_d = 2.0 * (1.0 - np.sqrt(1.0 - a/3.0))
        # calculate tessalation of mixed edges
        for i in range(31):
            j1 = 12
            j2 = 32
            a = ai_d
            if i > 12:
                j1 = i+1
                a = adod
            for j in range(j1, j2):
                d = sum((vertices[i] - vertices[j])**2)
                if abs(a-d) > DP_TOL:
                    continue
                for tl in range(1, tess):
                    vertices.append(_divarc(vertices[i], vertices[j], tl, tess))

        # calculate tessalation of pentakisdodecahedron faces
        for i in range(12):
            for j in range(12, 31):
                d = sum((vertices[i] - vertices[j])**2)
                if abs(ai_d-d) > DP_TOL:
                    continue
                for k in range(j+1, 32):
                    d_ik = sum((vertices[i] - vertices[k])**2)
                    d_jk = sum((vertices[j] - vertices[k])**2)
                    if (abs(ai_d - d_ik) > DP_TOL) or (abs(adod - d_jk) > DP_TOL):
                        continue
                    for tl in range(1, tess-1):
                        ji = _divarc(vertices[j], vertices[i], tl, tess)
                        ki = _divarc(vertices[k], vertices[i], tl, tess)
                        for tl2 in range(1, tess-tl):
                            ij = _divarc(vertices[i], vertices[j], tl2, tess)
                            kj = _divarc(vertices[k], vertices[j], tl2, tess)
                            ik = _divarc(vertices[i], vertices[k], tess-tl-tl2, tess)
                            jk = _divarc(vertices[j], vertices[k], tess-tl-tl2, tess)

                            xyz1 = _divarc(ki, ji, tl2, tess-tl)
                            xyz2 = _divarc(kj, ij, tl, tess-tl2)
                            xyz3 = _divarc(jk, ik, tl, tl+tl2)

                            x = xyz1 + xyz2 + xyz3
                            d = np.sqrt(sum(x**2))
                            vertices.append(x/d)

    return np.array(vertices)


def dotsphere(density):
    """Create a dot distribution over the unit sphere, choosing the most
    appropriate implementation based on the number of dots you request.
    
    Parameters
    ----------
    density : int
        Required number of dots on the unit sphere

    Returns
    -------
    dots : np.ndarray, shape=(N, 3), dtype=np.double
        Dots on the surface of the unit sphere. The number of dots will be
        at minimum equal to the `density` argument, but will be roughly two
        times larger.
    """
    
    i1 = 1
    i2 = 1
    while 10*i1*i1+2 < density:
        i1 += 1

    while 30*i2*i2+2 < density:
        i2 += 1

    if 10*i1*i1-2 < 30*i2*i2-2:
        return dotsphere1(density)
    else:
        return dotsphere2(density)


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


if __name__ == '__main__':
    for i in range(12, 80, 12):
        verts1 = dotsphere1(i)
        verts2 = dotsphere2(i)
