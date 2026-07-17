#include <iostream>
#include "ramses/core/Simulation.hpp"
#include "ramses/core/MpiManager.hpp"
#include "ramses/utils/Logger.hpp"

int main(int argc, char** argv) {
    auto& mpi = ramses::MpiManager::instance();
    mpi.init(argc, argv);
    
    ramses::Logger::init();

    RAMSES_INFO("\n{}", R"(
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
)");
    // Legacy read_params.f90:50 -- Working with nproc / ndim
    RAMSES_INFO(" Working with nproc = {:5d} for ndim = {:1d}", mpi.size(), NDIM);
    // Legacy read_params.f90:52-66 -- solver / nvar info
#ifdef MHD
    RAMSES_INFO(" Using solver = mhd with nvar = {:2d}", NDIM + 2);
#else
    RAMSES_INFO(" Using solver = hydro with nvar = {:2d}", NDIM + 2);
#endif
    RAMSES_INFO(" ");
    RAMSES_INFO(" compile date    = {}", RAMSES_BUILD_DATE);
    RAMSES_INFO(" compile command = {}", RAMSES_BUILD_COMMAND);
    RAMSES_INFO(" patch dir       = {}", RAMSES_PATCH_DIR);
    RAMSES_INFO(" remote repo     = {}", RAMSES_GIT_REPO);
    RAMSES_INFO(" local branch    = {}", RAMSES_GIT_BRANCH);
    RAMSES_INFO(" last commit     = {}", RAMSES_GIT_HASH);
    RAMSES_INFO(" ");

    if (argc < 2) {
        RAMSES_ERROR("Usage: {} <namelist.nml>", argv[0]);
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
