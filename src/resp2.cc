#include <stdio.h>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include "include/vdwsurface.h"
#define BOHR_TO_ANGSTROMS 0.52917721092

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace resp2 {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "RESP2"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_int("N_VDW_LAYERS", 4);
        options.add_double("VDW_SCALE_FACTOR", 1.4);
        options.add_double("VDW_INCREMENT", 0.2);
        options.add_double("VDW_POINT_DENSITY", 1.0);
        options.add_double("RESP_A", 0.005);
        options.add_double("RESP_B", 0.001);

    }

    return true;
}

extern "C" 
PsiReturnType resp2(Options& options)
{
    int print = options.get_int("PRINT");

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<BasisSet> basisset = wfn->basisset();
    boost::shared_ptr<Molecule> mol = basisset->molecule();
    boost::shared_ptr<IntegralFactory> integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset, basisset, basisset, basisset));
    boost::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));
    boost::shared_ptr<OEProp> oeprop = boost::shared_ptr<OEProp>(new OEProp());

    int n_atoms = mol->natom();    
    int nbf = basisset->nbf();

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
            nuc += mol->Z(atom1) / points[i].distance(mol->xyz(atom1));

        esp_values.push_back(nuc+elec);
        fprintf(outfile, "     %8.5f %8.5f %8.5f    %16.12f\n", points[i][0], points[i][1], points[i][2], nuc+elec);
    }
    fprintf(outfile, " ---------------------------------------------------\n\n");
    fflush(outfile);
    
    return Success;
}

}} // End namespaces

