#ifndef RAMSES_PARTICLE_SOLVER_HPP
#define RAMSES_PARTICLE_SOLVER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Handles particle dynamics (movement, force interpolation).
 */
class ParticleSolver {
public:
    ParticleSolver(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {}
    virtual ~ParticleSolver();

    /**
     * @brief Assigns particle mass to the grid density field (Cloud-In-Cell).
     */
    virtual void assign_mass(int ilevel);

    /**
     * @brief Pushes particles using current force field.
     */
    virtual void move_fine(int ilevel, real_t dt);

    /**
     * @brief Re-links particles to their current leaf cells.
     */
    virtual void relink();

    /**
     * @brief Assigns particle mass to the density grid (CIC).
     */
    virtual void assign_mass_fine(int ilevel);

    /**
     * @brief Performs MPI particle exchange between ranks.
     */
    virtual void exchange_particles();

    struct ParticlePacket {
        real_t xp[3];
        real_t vp[3];
        real_t mp;
        int idp;
        int levelp;
        uint8_t family;
        uint8_t tag;
    };

protected:
    virtual void move_particles(const std::vector<int>& ind_part, real_t dt);
    virtual int find_cell_by_coords(real_t x, real_t y, real_t z, int* level = nullptr);

    AmrGrid& grid_;
    Config& config_;
};

} // namespace ramses

#endif // RAMSES_PARTICLE_SOLVER_HPP
