#include <map>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/respfit.h"

using namespace boost::numeric;

//**************************************************************************//
// Resp fitting optimization
//**************************************************************************//

double resp_objective(const std::vector<double> &std_charges, std::vector<double> &grad, void *my_func_data) {
    Optdata* data = static_cast<Optdata*>(my_func_data);
    double a = data->resp_a;
    double b = data->resp_b;

    ublas::vector<double> charges(std_charges.size());
    ublas::vector<double> esp_values(data->esp_values.size());
    std::copy(std_charges.begin(), std_charges.end(), charges.begin());
    std::copy(data->esp_values.begin(), data->esp_values.end(), esp_values.begin());

    for (int i=0; i < std_charges.size(); i++)
        if (boost::math::isinf(std_charges[i]) || boost::math::isnan(std_charges[i]))
            printf("BADDDDD");

    // predicted esp values at the grid points based on the point charge model
    ublas::vector<double> esp_error = esp_values - ublas::prod(data->invr, charges);

    // figure of merit for how well the predicted charges match the actual
    double chi2_esp = ublas::norm_1(ublas::element_prod(esp_error, esp_error));

    // hyperbolic restraint term  a*sum(sqrt(q**2 + b**2)-b)
    double chi2_rstr = 0;
    for (size_t i = 0; i < charges.size(); i++)
        chi2_rstr += a * sqrt(charges[i]*charges[i] + b*b) - b;

    printf("Objective %f\n", chi2_esp + chi2_rstr);
    return chi2_esp + chi2_rstr;
}

double resp_constraint(const std::vector<double> &charges, std::vector<double> &grad, void* my_func_data) {
    // constraint function that should be equal to zero
    Optdata* data = static_cast<Optdata*>(my_func_data);

    // enforce that the total charge should be zero
    double total_charge = 0;
    for (size_t i = 0; i < charges.size(); i++)
        total_charge += charges[i];

    // we want error to be zero
    double error = total_charge*total_charge;

    // put the charge groups into a map
    std::map<int, std::vector<int> > groups;
    for (int i = 0; i < data->charge_groups.size(); i++) {
        int item = data->charge_groups[i];
        if (groups.find(item) == groups.end())
            groups[item] = std::vector<int>();
        groups[item].push_back(i);
    }

    for (std::map<int, std::vector<int> >::iterator it = groups.begin(); it != groups.end(); it++) {
        std::vector<int> items = it->second;
        for (int i = 0; i < items.size(); i++) {
            double diff = charges[items[i]] - charges[items[0]];
            error += diff*diff;
        }
    }
    return error;
}

