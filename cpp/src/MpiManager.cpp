#include "ramses/MpiManager.hpp"
#include <iostream>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

void MpiManager::init(int argc, char** argv) {
#ifdef RAMSES_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
#else
    rank_ = 0;
    size_ = 1;
#endif
    
    if (is_master()) {
        std::cout << "[MPI] Initialized. Rank=" << rank_ << " Size=" << size_ << std::endl;
    }
}

void MpiManager::finalize() {
#ifdef RAMSES_USE_MPI
    MPI_Finalize();
#endif
}

void MpiManager::barrier() {
#ifdef RAMSES_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

} // namespace ramses
