"""
Calculation of the fused-sphere van der waals surface of a molecule.
"""

#############################################################################
# Imports
#############################################################################

from __future__ import division
import itertools
import numpy as np

from dotsphere import dotsphere, plot_vertices

#############################################################################
# Globals
#############################################################################

BONDI_RADII = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.8,
    'S': 1.8,
    'Cl': 1.75}

#############################################################################
# Functions
#############################################################################

def vdw_surface(coordinates, elements, scale_factor=1.4, density=1.0):
    """Compute points on the VDW surface of a molecule
    
    Parameters
    ----------
    coordinates : np.ndarray, shape=(n_atoms, 3)
        The cartesian coordinates of the nuclei, in units of ANGSTROMS
    elements : list, shape=(n_atoms)
        The element symbols (C, H, etc) for each of the nuceli
    scale_factor : float
        The points on the molecular surface are set at a distance of
        scale_factor * vdw_radius away from each of the atoms.
    density : float
        The (approximate) number of points to generate per square angstrom
        of surface area. 1.0 is the default recommended by Kollman & Singh.
    """
    n_atoms = len(coordinates)
    if len(coordinates) != len(elements):
        raise ValueError('The number of atoms in coordinates (%d) and elements'
                         '(%d) do not match' % (len(coordinates), len(elements)))
    
    try:
        radii = [BONDI_RADII[e]*scale_factor for e in elements]
    except KeyError:
        raise KeyError('element must be one of %s' % str(BONDI_RADII.keys()))
    
    surfacepoints = []
    compute = 0
    
    for i in range(n_atoms):
        # this could be optimized -- we don't need to compute the dotsphere.
        # at maximum we only need to compute it once per unique element / radius
        unit_dots = dotsphere(density * ((4/3) * np.pi * radii[i]**3))
        dots = coordinates[i] + radii[i] * unit_dots
    
        # all of the atoms that i is close to
        neighbors = []
        neighbordistance = []
        
        for j in range(n_atoms):
            if j == i:
                continue
            d = sum((coordinates[i] - coordinates[j])**2)
            if d < (radii[i] + radii[j])**2:
                neighbors.append(j)
                neighbordistance.append(d)
        
        # sort the neighbors by distance to j. this speeds up finding
        # inaccessible dots and rejecting them quickly
        neighbors = np.array(neighbors)[np.argsort(neighbordistance)]
        
        for k in range(len(dots)):
            accessible = True
            
            for j in neighbors:
                # if the dot is within the vdw radius of the atom j
                compute += 1
                if sum((dots[k] - coordinates[j])**2) < radii[j]**2:
                    jstart = j
                    accessible = False
                    break
            
            if accessible:
                surfacepoints.append(dots[k])

    return np.array(surfacepoints)


if __name__ == '__main__':
    # plot the surface for ala2, just to make sure it looks good
    
    import mdtraj as md
    t = md.load('ala2.pdb')
    elements = [a.element.symbol for a in t.topology.atoms]
    
    # convert nm to A
    xyz = t.xyz[0]*10

    plot_vertices(vdw_surface(xyz, elements, density=2))
    