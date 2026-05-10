#include "ramses/Parameters.hpp"

namespace ramses {
namespace params {

    bool verbose = false;
    bool hydro = false;
    bool pic = false;
    bool poisson = false;
    bool cosmo = false;

    int nx = 1, ny = 1, nz = 1;
    int levelmin = 1;
    int nlevelmax = 1;
    int ngridmax = 0;
    real_t boxlen = 1.0;
    int iriemann = 1; // Default to LLF

    int nstepmax = 1000000;
    real_t trestart = 0.0;

    real_t units_length = 1.0;
    real_t units_density = 1.0;
    real_t units_time = 1.0;
    real_t units_velocity = 1.0;
    real_t units_mass = 1.0;
    real_t units_energy = 1.0;
    real_t units_pressure = 1.0;
    real_t units_number_density = 1.0;

    bool barotropic_eos = false;
    std::string barotropic_eos_form = "isothermal";
    real_t polytrope_rho = 1e-15;
    real_t polytrope_index = 1.4;
    real_t T_eos = 10.0;
    real_t mu_gas = 1.0;

} // namespace params
} // namespace ramses
