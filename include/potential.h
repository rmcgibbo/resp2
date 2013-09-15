#ifndef _RESP_POTENTIAL_H
#define _RESP_POTENTIAL_H

#include <stdio.h>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <nlopt.hpp>

std::vector<double> calculate_esp_at_points(std::vector<Vector3> points);

#endif