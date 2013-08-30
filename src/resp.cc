#include <stdio.h>

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace resp {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "RESP"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_int("N_SHELLS", 4);
    }

    return true;
}

extern "C" 
PsiReturnType resp(Options& options)
{
    int print = options.get_int("PRINT");
    
    printf("n_shells: %d", options.get_int("N_SHELLS"));
    /* Your code goes here */
    printf("resp.cc: resp\n");
    return Success;
}

}} // End namespaces

