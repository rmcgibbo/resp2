#ifndef _dotsphere_h
#define _dotsphere_h

#include <vector>
// this is also included in the psi4 code as libmints/vector3.h
#include "vector3.h"

using namespace std;
using namespace psi;

vector<Vector3> dotsphere(int density);

#endif
