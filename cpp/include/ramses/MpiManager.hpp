#ifndef RAMSES_MPI_MANAGER_HPP
#define RAMSES_MPI_MANAGER_HPP

#include <string>

namespace ramses {

/**
 * @brief Singleton-like manager for MPI communication.
 * 
 * Replicates amr/mpi_mod.f90 functionality.
 */
class MpiManager {
public:
    static MpiManager& instance() {
        static MpiManager inst;
        return inst;
    }

    void init(int argc, char** argv);
    void finalize();

    int rank() const { return rank_; }
    int size() const { return size_; }
    bool is_master() const { return rank_ == 0; }

    void barrier();

private:
    MpiManager() : rank_(0), size_(1) {}
    int rank_;
    int size_;
};

} // namespace ramses

#endif // RAMSES_MPI_MANAGER_HPP
