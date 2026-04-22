#ifndef RAMSES_PARTICLE_SOLVER_HPP
#define RAMSES_PARTICLE_SOLVER_HPP

#include "AmrGrid.hpp"
#include "ParticleSystem.hpp"

namespace ramses {

/**
 * @brief Handles particle dynamics (movement, force interpolation).
 */
class ParticleSolver {
public:
    ParticleSolver(AmrGrid& grid, ParticleSystem& ps) : grid_(grid), ps_(ps) {}

    /**
     * @brief Assigns particle mass to the grid density field (Cloud-In-Cell).
     */
    void assign_mass(int ilevel);

    /**
     * @brief Moves particles at a specific level.
     */
    void move_fine(int ilevel, real_t dt);

private:
    /**
     * @brief Core advection logic for a batch of particles.
     */
    void move_particles(const std::vector<int>& ind_part, real_t dt);

    AmrGrid& grid_;
    ParticleSystem& ps_;
};

} // namespace ramses

#endif // RAMSES_PARTICLE_SOLVER_HPP
