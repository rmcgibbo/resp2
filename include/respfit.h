#ifndef _RESP_RESPFIT_H
#define _RESP_RESPFIT_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

typedef struct {
    ublas::matrix<double> invr;
    std::vector<double> esp_values;
    std::vector<int> charge_groups;
    double resp_a;
    double resp_b;
    int n_iterations;
    std::vector<double> mulliken_charges;
} Optdata;

/*
 * resp_objective computes 
 */
double resp_objective(const std::vector<double> &std_charges, std::vector<double> &grad, void *my_func_data);

/*
 * The constraints
 */
double resp_constraint(const std::vector<double> &charges, std::vector<double> &grad, void* my_func_data);

#endif