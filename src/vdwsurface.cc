//  Copyright (c) 2013, Robert McGibbon
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////
// Calculation of the fused-sphere van der waals surface of a molecule.
//////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <set>

#include "include/dotsphere.h"
// this is also included in the psi4 code as libmints/vector3.h
#include "include/vector3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace psi;


// A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001
static map<string, double> create_bondi_radii() {
    map<string, double> m;
    m["H"] = 1.2;
    m["C"] = 1.7;
    m["N"] = 1.55;
    m["O"] = 1.52;
    m["F"] = 1.47;
    m["P"] = 1.8;
    m["S"] = 1.8;
    m["Cl"] = 1.75;

    m["Ar"] = 1.88;
    m["As"] = 1.85;
    m["Br"] = 1.85;
    m["Cd"] = 1.62;
    m["Cu"] = 1.4;
    m["Ga"] = 1.87;
    m["Au"] = 1.66;
    m["He"] = 1.4;
    m["In"] = 1.93;
    m["I"] = 1.98;
    m["Kr"] = 2.02;
    m["Pb"] = 2.02;
    m["Li"] = 1.82;
    m["Mg"] = 1.73;
    m["Hg"] = 1.70;
    m["Ne"] = 1.54;
    m["Ni"] = 1.64;
    m["Pd"] = 1.63;
    m["Pt"] = 1.8;
    m["K"] = 2.75;
    m["Se"] = 1.90;
    m["Si"] = 2.1;
    m["Ag"] = 1.9;
    m["Na"] = 2.27;
    m["Te"] = 2.06;
    m["Tl"] = 1.96;
    m["Sn"] = 2.17;
    m["U"] = 1.86;
    m["Xe"] = 2.16;
    m["Zn"] = 1.37;
    return m;
}
static map<string, double> BONDI_RADII = create_bondi_radii();


// #############################################################################
// # Functions
// #############################################################################
//
vector<Vector3> vdw_surface(vector<Vector3> coordinates, vector<string> elements,
                            double scale_factor, double density) {
    // Compute points on the VDW surface of a molecule
    //
    // Parameters
    // ----------
    // coordinates : np.ndarray, shape=(n_atoms, 3)
    //     The cartesian coordinates of the nuclei, in units of ANGSTROMS
    // elements : list, shape=(n_atoms)
    //     The element symbols (C, H, etc) for each of the nuceli
    // scale_factor : float
    //     The points on the molecular surface are set at a distance of
    //     scale_factor * vdw_radius away from each of the atoms.
    // density : float
    //     The (approximate) number of points to generate per square angstrom
    //     of surface area. 1.0 is the default recommended by Kollman & Singh.
    if (coordinates.size() != elements.size()) {
        fprintf(stderr, "coordinate.size doesnot match elements.size");
        exit(1);
    }

    vector<Vector3> surfacepoints;
    vector<double> radii;

    for (size_t i = 0; i < elements.size(); i++) {
        // todo: check for error if element not in BONDI_RADII table
        if (BONDI_RADII.find(elements[i]) != BONDI_RADII.end()) {
            radii.push_back(BONDI_RADII[elements[i]] * scale_factor);
        } else {
            fprintf(stderr, "%s is not a supported element", elements[i].c_str());
            exit(1);
        }
    }

    for (size_t i = 0; i < coordinates.size(); i++) {
        // this could be optimized -- we don't need to compute the dotsphere.
        // at maximum we only need to compute it once per unique element / radius
        vector<Vector3> dots = dotsphere(density * ((4.0/3.0) * M_PI * radii[i]*radii[i]*radii[i]));

        for (size_t j = 0; j < dots.size(); j++)
            dots[j] = coordinates[i] + radii[i]*dots[j];

        // all of the atoms that i is close to
        typedef std::pair<double, size_t> Neighbor;
        vector<Neighbor> neighbors;
        for (size_t j = 0; j < coordinates.size(); j++) {
            if (i == j)
                continue;
            double d = (coordinates[i] - coordinates[j]).norm();
            if (d < (radii[i] + radii[j])) {
                neighbors.push_back(make_pair(d, j));
            }
        }
        sort(neighbors.begin(), neighbors.end());

        for (size_t k = 0; k < dots.size(); k++) {
            int accessible = 1;
            for (vector<Neighbor>::iterator it = neighbors.begin() ; it != neighbors.end(); ++it) {
                if ((coordinates[(*it).second] - dots[k]).norm() < radii[(*it).second]) {
                    accessible = 0;
                    break;
                }
            }
            if (accessible)
                surfacepoints.push_back(dots[k]);
        }
    }

    return surfacepoints;
}

// extern "C" void vdw_surface(double* coordinates, char* elements, int n_elements,
//                             double scale_factor, double density, double* out, int* n_out) {
//     vector<string> elements_;
//     vector<Vector3> coordinates_;
// 
//     istringstream iss(elements);
//     string s;
//     for (int i = 0; i < n_elements; i++) {
//         getline(iss, s, ' ');
//         elements_.push_back(string(s));
//         coordinates_.push_back(Vector3(coordinates[3*i+0], coordinates[3*i+1], coordinates[3*i+2]));
//     }
// 
//     vector<Vector3> points = vdw_surface(coordinates_, elements_, scale_factor, density);
//     for (size_t i = 0; i < points.size(); i++) {
//         out[3*i+0] = points[i][0];
//         out[3*i+1] = points[i][1];
//         out[3*i+2] = points[i][2];
//     }
//     *n_out = points.size();
// }
