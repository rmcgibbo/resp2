#ifndef _vdwsurface_h
#define _vdwsurface_h

#include <vector>
#include "libmints/vector3.h"

using namespace std;
using namespace psi;

vector<Vector3> vdw_surface(vector<Vector3> coordinates, vector<string> elements,
  double scale_factor, double density);;
                              
#endif
