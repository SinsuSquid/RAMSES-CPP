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
    extern int iriemann; // 1=llf, 2=hllc

    // Time step parameters
    extern int nstepmax;
    extern real_t trestart;

    // Physical units (cgs)
    extern real_t units_length;
    extern real_t units_density;
    extern real_t units_time;
    extern real_t units_velocity;
    extern real_t units_mass;
    extern real_t units_energy;
    extern real_t units_pressure;
    extern real_t units_number_density;

    // Barotropic EOS parameters
    extern bool barotropic_eos;
    extern std::string barotropic_eos_form;
    extern real_t polytrope_rho;
    extern real_t polytrope_index;
    extern real_t T_eos;
    extern real_t mu_gas;

    // Constants derived from NDIM
    constexpr int twotondim = (1 << NDIM);
    constexpr int threetondim = 1; // Needs to be calculated based on NDIM
    constexpr int twondim = 2 * NDIM;

} // namespace params
} // namespace ramses

#endif // RAMSES_PARAMETERS_HPP
