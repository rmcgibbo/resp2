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
#include "include/potential.h"
#include "include/respfit.h"

#define BOHR_TO_ANGSTROMS 0.52917721092
#define CONSTRAINT_TOLERANCE 1e-8
#define CHARGE_TOLERANCE 1e-8

INIT_PLUGIN

using namespace boost;
using namespace boost::numeric;

namespace psi{ namespace resp2 {

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
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if (!wfn) {
        fprintf(outfile, "\n ---------------------------------------------------\n");
        fprintf(outfile, " You must first compute the wavefunction (scf) before\n");
        fprintf(outfile, " running the resp plugin\n");
        fprintf(outfile, " ---------------------------------------------------\n");
            return Failure;
    }

    boost::shared_ptr<Molecule> mol = wfn->basisset()->molecule();
    int n_atoms = mol->natom();

    // Make sure that the charge groups are reasonable
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

    // Get the coordinates of the nuclei in Angstroms
    double point_density = options.get_double("VDW_POINT_DENSITY");
    std::vector<std::string> symbols;
    std::vector<Vector3> coordinates;
    for (int i = 0; i < n_atoms; i++) {
        symbols.push_back(mol->symbol(i));
        // the vdwsurface code expects coordinates in angstroms
        coordinates.push_back(Vector3(mol->xyz(i)) * BOHR_TO_ANGSTROMS);
    }

    // Get the points at which we're going to calculate the ESP surface
    std::vector<Vector3> points;
    for (int i = 0; i < options.get_int("N_VDW_LAYERS"); i++) {
        double scale_factor = options.get_double("VDW_SCALE_FACTOR") + i * options.get_double("VDW_INCREMENT");
        std::vector<Vector3> this_shell = vdw_surface(coordinates, symbols, scale_factor, point_density);
        points.insert(points.end(), this_shell.begin(), this_shell.end());
    }

    // Convert these ESP points back into BOHRS
    for (size_t i = 0; i < points.size(); i++)
        points[i] = points[i] / BOHR_TO_ANGSTROMS;
    // And then calculate the ESP values there
    std::vector<double> esp_values = calculate_esp_at_points(points);

    // Build a matrix of the inverse distance from each ESP point to each nucleus
    ublas::matrix<double> invr (points.size(), n_atoms);
    for (size_t i = 0; i < invr.size1(); i++)
        for (size_t j = 0; j < invr.size2(); j++)
            invr(i, j) = 1.0 / points[i].distance(mol->xyz(j));
        
    // Run the optimization using NLOPT
    nlopt::opt opt(nlopt::LN_COBYLA, n_atoms);
    Optdata data;
    data.invr = invr;
    data.esp_values = esp_values;
    data.resp_a = options.get_double("RESP_A");
    data.resp_b = options.get_double("RESP_B");
    data.charge_groups = charge_groups;
    data.n_iterations = 0;
    
    opt.set_min_objective(resp_objective, &data);
    opt.set_xtol_abs(CHARGE_TOLERANCE);
    opt.add_equality_constraint(resp_constraint, &data, CONSTRAINT_TOLERANCE);
    std::vector<double> charges(n_atoms, 0.0);

    double opt_f = 0;
    nlopt::result result = opt.optimize(charges, opt_f);
    
    // Print the parameters to disk
    fprintf(outfile, "\n ---------------------------------------------------\n");
    fprintf(outfile, " RESTRAINED ELECTROSTATIC POTENTIAL PARAMETERS\n");
    fprintf(outfile, " ---------------------------------------------------\n");
    fprintf(outfile, " N_VDW_LAYERS:       %d\n", options.get_int("N_VDW_LAYERS"));
    fprintf(outfile, " VDW_SCALE_FACTOR:   %.3f\n", options.get_double("VDW_SCALE_FACTOR"));
    fprintf(outfile, " VDW_INCREMENT:      %.3f\n", options.get_double("VDW_INCREMENT"));
    fprintf(outfile, " VDW_POINT_DENSITY:  %.3f\n", options.get_double("VDW_POINT_DENSITY"));
    fprintf(outfile, " RESP_A:             %.4f\n", options.get_double("RESP_A"));
    fprintf(outfile, " RESP_B:             %.4f\n", options.get_double("RESP_B"));
    fprintf(outfile, " CHARGE_GROUPS:      [");
    for (size_t i = 0; i < charge_groups.size()-1; i++)
        fprintf(outfile, "%d, ", charge_groups[i]);
    fprintf(outfile, "%d]\n", charge_groups[charge_groups.size()-1]);
    fprintf(outfile, " ---------------------------------------------------\n");
    
    // Print the results to disk
    fprintf(outfile, "\n ----------------------------------------------\n");
    fprintf(outfile, " RESTRAINED ELECTROSTATIC POTENTIAL CHARGES\n");
    fprintf(outfile, "   Center  Symbol  RESP Charge (a.u.)\n");
    fprintf(outfile, " ----------------------------------------------\n");
    for (int i = 0; i < charges.size(); i++)
        fprintf(outfile, "   %5d    %2s     %8.5f\n", i+1, symbols[i].c_str(), charges[i]);
    fprintf(outfile, " ----------------------------------------------\n");

    return Success;
}

}} // End namespaces
