#include <stdio.h>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

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
    
    printf("n_shells: %d\n", options.get_int("N_VDW_LAYERS"));
    /* Your code goes here */
    printf("resp.cc2: resp2()\n");
    return Success;
}

}} // End namespaces

