#include <iostream>
#include "ramses/core/Simulation.hpp"
#include "ramses/core/MpiManager.hpp"

int main(int argc, char** argv) {
    auto& mpi = ramses::MpiManager::instance();
    mpi.init(argc, argv);

    if (mpi.is_master()) {
        std::cout << R"(
 ██████╗  █████╗ ███╗   ███╗███████╗███████╗███████╗
 ██╔══██╗██╔══██╗████╗ ████║██╔════╝██╔════╝██╔════╝
 ██████╔╝███████║██╔████╔██║███████╗█████╗  ███████╗
 ██╔══██╗██╔══██║██║╚██╔╝██║╚════██║██╔══╝  ╚════██║
 ██║  ██║██║  ██║██║ ╚═╝ ██║███████║███████╗███████║
 ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝╚══════╝╚══════╝

                ::  ::       CCCCCC   PPPPPPPP  PPPPPPPP
                ::  ::      CC    CC  PP    PP  PP    PP
                            CC        PPPPPPPP  PPPPPPPP
                ::  ::      CC    CC  PP        PP
                ::  ::       CCCCCC   PP        PP      ;
)" << std::endl;
        // Legacy read_params.f90:50 -- Working with nproc / ndim
        printf(" Working with nproc = %5d for ndim = %1d\n", mpi.size(), NDIM);
        // Legacy read_params.f90:52-66 -- solver / nvar info
#ifdef MHD
        printf(" Using solver = mhd with nvar = %2d\n", NDIM + 2);
#else
        printf(" Using solver = hydro with nvar = %2d\n", NDIM + 2);
#endif
        printf(" \n");
        printf(" compile date    = %s\n", RAMSES_BUILD_DATE);
        printf(" compile command = %s\n", RAMSES_BUILD_COMMAND);
        printf(" patch dir       = %s\n", RAMSES_PATCH_DIR);
        printf(" remote repo     = %s\n", RAMSES_GIT_REPO);
        printf(" local branch    = %s\n", RAMSES_GIT_BRANCH);
        printf(" last commit     = %s\n", RAMSES_GIT_HASH);
        printf(" \n");
    }

    if (argc < 2) {
        if (mpi.is_master()) {
            std::cout << "Usage: " << argv[0] << " <namelist.nml>" << std::endl;
        }
        mpi.finalize();
        return 1;
    }

    std::string nml_path = argv[1];
    
    ramses::Simulation sim;
    sim.initialize(nml_path);
    sim.run();

    mpi.finalize();
    return 0;
}
