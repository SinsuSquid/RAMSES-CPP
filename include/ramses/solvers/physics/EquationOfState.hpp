#ifndef RAMSES_EQUATION_OF_STATE_HPP
#define RAMSES_EQUATION_OF_STATE_HPP

#include "Types.hpp"
#include <vector>
#include <string>

namespace ramses {

/**
 * @brief Unified helper class for Equation of State (EoS) and thermodynamic calculations.
 * 
 * Centralizes the calculation of pressures, sound speeds, and enthalpies
 * across hydro, MHD, and relativistic solvers.
 */
class EquationOfState {
public:
    /**
     * @brief Compute thermal pressure for ideal gas EoS.
     */
    static real_t get_pressure(real_t density, real_t thermal_energy_density, real_t gamma);

    /**
     * @brief Compute pressure for barotropic EoS from density.
     */
    static real_t get_barotropic_pressure(real_t density);

    /**
     * @brief Compute sound speed squared (cs^2) for hydro/MHD.
     */
    static real_t get_cs2(real_t density, real_t pressure, real_t gamma,
                          const real_t q[] = nullptr, int nener = 0,
                          const std::vector<real_t>& gamma_rad = {});

    /**
     * @brief Compute specific enthalpy (h) for relativistic hydrodynamics (RHD).
     */
    static real_t get_rhd_enthalpy(real_t density, real_t pressure, real_t gamma, const std::string& eos_type);

    /**
     * @brief Compute sound speed (cs) for RHD.
     */
    static real_t get_rhd_sound_speed(real_t density, real_t pressure, real_t gamma, const std::string& eos_type);

    /**
     * @brief Compute pressure for RHD when momentum M == 0.
     */
    static real_t get_rhd_pressure_zero_momentum(real_t D, real_t E, real_t gamma, const std::string& eos_type);

    /**
     * @brief Compute pressure for RHD from the auxiliary variable xsi and density rho.
     */
    static real_t get_rhd_pressure_from_xsi(real_t xsi, real_t rho, real_t gamma, const std::string& eos_type);
};

} // namespace ramses

#endif // RAMSES_EQUATION_OF_STATE_HPP
