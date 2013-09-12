#include <stdio.h>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <nlopt.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "include/vdwsurface.h"
#define BOHR_TO_ANGSTROMS 0.52917721092
#define CONSTRAINT_TOLERANCE 1e-8
#define CHARGE_TOLERANCE 1e-8

INIT_PLUGIN

using namespace boost;
using namespace boost::numeric;

namespace psi{ namespace resp2 {

//**************************************************************************//
// Resp fitting optimization
//**************************************************************************//


typedef struct {
    ublas::matrix<double> invr;
    std::vector<double> esp_values;
    std::vector<int> charge_groups;
    double resp_a;
    double resp_b;
} Optdata;


double resp_objective(const std::vector<double> &std_charges, std::vector<double> &grad, void *my_func_data) {
    Optdata* data = static_cast<Optdata*>(my_func_data);
    double a = data->resp_a;
    double b = data->resp_b;

    ublas::vector<double> charges(std_charges.size());
    ublas::vector<double> esp_values(data->esp_values.size());
    std::copy(std_charges.begin(), std_charges.end(), charges.begin());
    std::copy(data->esp_values.begin(), data->esp_values.end(), esp_values.begin());

    // predicted esp values at the grid points based on the point charge model
    ublas::vector<double> esp_error = esp_values - ublas::prod(data->invr, charges);

    // figure of merit for how well the predicted charges match the actual
    double chi2_esp = ublas::norm_1(ublas::element_prod(esp_error, esp_error));

    // hyperbolic restraint term  a*sum(sqrt(q**2 + b**2)-b)
    double chi2_rstr = 0;
    for (size_t i = 0; i < charges.size(); i++)
        chi2_rstr += a * sqrt(charges[i]*charges[i] + b*b) - b;

    //printf("Objective %f\n", chi2_esp + chi2_rstr);
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


//**************************************************************************//
// Interface with psi4
//**************************************************************************//


extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "RESP2"|| options.read_globals()) {
        options.add_int("N_VDW_LAYERS", 4);
        options.add_double("VDW_SCALE_FACTOR", 1.4);
        options.add_double("VDW_INCREMENT", 0.2);
        options.add_double("VDW_POINT_DENSITY", 1.0);
        options.add_double("RESP_A", 0.005);
        options.add_double("RESP_B", 0.001);
        options.add_array("CHARGE_GROUPS");

    }

    return true;
}

//**************************************************************************//
// "main" function for the plugin
//**************************************************************************//


extern "C"
PsiReturnType resp2(Options& options) {
    fprintf(outfile, "\n ---------------------------------------------------\n");
    fprintf(outfile, " RESTRAINED ELECTROSTATIC POTENTIAL CHARGE FITTING\n");
    fprintf(outfile, " ---------------------------------------------------\n");

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<BasisSet> basisset = wfn->basisset();
    boost::shared_ptr<Molecule> mol = basisset->molecule();
    boost::shared_ptr<IntegralFactory> integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset, basisset, basisset, basisset));
    boost::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));
    boost::shared_ptr<OEProp> oeprop = boost::shared_ptr<OEProp>(new OEProp());

    int n_atoms = mol->natom();
    int nbf = basisset->nbf();

    std::vector<int> charge_groups = options.get_int_vector("CHARGE_GROUPS");
    if (charge_groups.size() == 0) {
        fprintf(outfile, "RESP: Using default charge groups. Every atom will be assigned a unique charge\n");
        for (int i = 0; i < n_atoms; i++)
            charge_groups.push_back(i);
    } else if (charge_groups.size() != n_atoms) {
        fprintf(outfile, "RESP: FATAL ERROR\n");
        fprintf(outfile, "CHARGE_GROUPS must be a list of integers of size equal to the number of atoms\n");
        exit(1);
    }

    double point_density = options.get_double("VDW_POINT_DENSITY");
    std::vector<std::string> symbols;
    std::vector<Vector3> coordinates;
    for (int i = 0; i < n_atoms; i++) {
        symbols.push_back(mol->symbol(i));
        // the vdwsurface code expects coordinates in angstroms
        coordinates.push_back(Vector3(mol->xyz(i)) * BOHR_TO_ANGSTROMS);
    }

    // the the points at which we're going to calculate the ESP surface
    std::vector<Vector3> points;
    std::vector<double> esp_values;
    for (int i = 0; i < options.get_int("N_VDW_LAYERS"); i++) {
        double scale_factor = options.get_double("VDW_SCALE_FACTOR") + i * options.get_double("VDW_INCREMENT");
        std::vector<Vector3> this_shell = vdw_surface(coordinates, symbols, scale_factor, point_density);
        points.insert(points.end(), this_shell.begin(), this_shell.end());
    }

    // convert the points back to bohrs
    for (size_t i = 0; i < points.size(); i++)
        points[i] = points[i] / BOHR_TO_ANGSTROMS;

    SharedMatrix Dtot = oeprop->Da_ao();
    if (wfn->same_a_b_dens()) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(oeprop->Db_ao());
    }


    // compute the electrostatic potential at each of the points
    // this code probably should be in a different function?
    fprintf(outfile, "\n Electrostatic potentials at van der Waals shells:\n");
    fprintf(outfile, " ---------------------------------------------------\n");
    fprintf(outfile, "   x   y   z  Electrostatic Potential (a.u.)\n");
    fprintf(outfile, " ---------------------------------------------------\n");
    for (size_t i = 0; i < points.size(); i++) {
        std::stringstream s;
        SharedMatrix ints(new Matrix(s.str(), nbf, nbf));
        epot->compute(ints, points[i]);
        double elec = Dtot->vector_dot(ints);
        double nuc = 0.0;

        for (int atom1 = 0; atom1 < n_atoms; atom1++)
            nuc += (mol->Z(atom1) / points[i].distance(mol->xyz(atom1)));

        esp_values.push_back(nuc+elec);
        fprintf(outfile, "     %8.5f %8.5f %8.5f    %16.12f\n", points[i][0], points[i][1], points[i][2], nuc+elec);
    }
    fprintf(outfile, " ---------------------------------------------------\n\n");
    fflush(outfile);

    printf("Fitting ESP...\n");

    ublas::matrix<double> invr (points.size(), n_atoms);
    for (size_t i = 0; i < invr.size1(); i++)
        for (size_t j = 0; j < invr.size2(); j++)
            invr(i, j) = 1.0 / points[i].distance(mol->xyz(j));

    nlopt::opt opt(nlopt::LN_COBYLA, n_atoms);
    Optdata data;
    data.invr = invr;
    data.esp_values = esp_values;
    data.resp_a = options.get_double("RESP_A");
    data.resp_b = options.get_double("RESP_B");
    data.charge_groups = charge_groups;

    opt.set_min_objective(resp_objective, &data);
    opt.set_xtol_abs(CHARGE_TOLERANCE);
    opt.add_equality_constraint(resp_constraint, &data, CONSTRAINT_TOLERANCE);

    std::vector<double> charges(n_atoms);
    for (size_t i = 0; i < n_atoms; i++)
        charges[i] = 0.0;

    double opt_f = 0;
    nlopt::result result = opt.optimize(charges, opt_f);

    for (int i = 0; i < charges.size(); i++)
        printf("%d: %f    ", i, charges[i]);
    printf("\n");

    return Success;
}



}} // End namespaces
