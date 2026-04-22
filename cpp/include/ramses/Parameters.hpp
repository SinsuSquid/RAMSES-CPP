#ifndef RAMSES_PARAMETERS_HPP
#define RAMSES_PARAMETERS_HPP

#include "Types.hpp"
#include <string>
#include <vector>

namespace ramses {
namespace params {

    // Run control flags
    extern bool verbose;
    extern bool hydro;
    extern bool pic;
    extern bool poisson;
    extern bool cosmo;

    // Mesh parameters
    extern int nx, ny, nz;
    extern int levelmin;
    extern int nlevelmax;
    extern int ngridmax;
    extern real_t boxlen;

    // Time step parameters
    extern int nstepmax;
    extern real_t trestart;

    // Constants derived from NDIM
    constexpr int twotondim = (1 << NDIM);
    constexpr int threetondim = 1; // Needs to be calculated based on NDIM
    constexpr int twondim = 2 * NDIM;

} // namespace params
} // namespace ramses

#endif // RAMSES_PARAMETERS_HPP
