#ifndef RAMSES_PARTICLE_SYSTEM_HPP
#define RAMSES_PARTICLE_SYSTEM_HPP

#include "Types.hpp"
#include <vector>
#include <string>

namespace ramses {

/**
 * @brief Manages N-body particles (Dark Matter, Stars, etc.).
 * 
 * Replicates arrays in pm/pm_commons.f90.
 */
class ParticleSystem {
public:
    ParticleSystem() = default;

    /**
     * @brief Allocates memory for a maximum number of particles.
     */
    void allocate(int npartmax, int ngridmax);

    // Particle data (AoB - Array of Builders / Struct of Arrays)
    std::vector<real_t> xp;       // Positions [npartmax * NDIM]
    std::vector<real_t> vp;       // Velocities [npartmax * NDIM]
    std::vector<real_t> mp;       // Masses [npartmax]
    std::vector<int> levelp;      // Level of particle [npartmax]
    std::vector<i8b_t> idp;       // Identity [npartmax]
    
    // Linked list for particles in cells
    std::vector<int> nextp;       // Next particle in list [npartmax]
    std::vector<int> prevp;       // Previous particle in list [npartmax]

    // Grid pointers
    std::vector<int> headp;       // Head particle in grid [ngridmax * twotondim]
    std::vector<int> tailp;       // Tail particle in grid [ngridmax * twotondim]
    std::vector<int> numbp;       // Number of particles in grid [ngridmax * twotondim]

    // Free list pointers
    int headp_free, tailp_free, numbp_free;

    int npartmax;
    int ngridmax;

    // 1-based utility accessors
    inline real_t& get_xp(int ipart, int idim) { return xp[(idim - 1) * npartmax + (ipart - 1)]; }
    inline real_t& get_vp(int ipart, int idim) { return vp[(idim - 1) * npartmax + (ipart - 1)]; }
};

} // namespace ramses

#endif // RAMSES_PARTICLE_SYSTEM_HPP
