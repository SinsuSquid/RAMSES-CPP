#include "ramses/solvers/physics/EquationOfState.hpp"
#include "ramses/core/Parameters.hpp"
#include "ramses/core/Constants.hpp"
#include <cmath>
#include <algorithm>

namespace ramses {

real_t EquationOfState::get_pressure(real_t density, real_t thermal_energy_density, real_t gamma) {
    real_t p = std::max(thermal_energy_density * (gamma - 1.0), density * 1e-10);
    // Cap pressure to prevent sound speed explosion in vacuum
    return std::min(p, density * 1e6);
}

real_t EquationOfState::get_barotropic_pressure(real_t density) {
    real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
    real_t temp_mu = T2;
    if (params::barotropic_eos_form == "polytrope") {
        temp_mu = T2 * std::pow(density * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
    } else if (params::barotropic_eos_form == "double_polytrope") {
        temp_mu = T2 * (1.0 + std::pow(density * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0));
    }
    return density * temp_mu / (params::units_velocity * params::units_velocity);
}

real_t EquationOfState::get_cs2(real_t density, real_t pressure, real_t gamma,
                                const real_t q[], int nener,
                                const std::vector<real_t>& gamma_rad) {
    if (params::barotropic_eos) {
        real_t T2 = params::T_eos / params::mu_gas * (constants::kB / constants::mH);
        real_t v2_unit = params::units_velocity * params::units_velocity;
        if (params::barotropic_eos_form == "isothermal") {
            return T2 / v2_unit;
        }
        if (params::barotropic_eos_form == "polytrope") {
            return params::polytrope_index * T2 * std::pow(density * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0) / v2_unit;
        }
        if (params::barotropic_eos_form == "double_polytrope") {
            real_t fac = std::pow(density * params::units_density / params::polytrope_rho, params::polytrope_index - 1.0);
            return T2 * (1.0 + params::polytrope_index * fac) / v2_unit;
        }
    }
    
    real_t cs2 = gamma * pressure;
    if (q != nullptr) {
        for (int ie = 0; ie < nener; ++ie) {
            cs2 += gamma_rad[ie] * q[NDIM + 2 + ie];
        }
    }
    return std::clamp(cs2 / std::max(density, 1e-10), (real_t)1e-2, (real_t)1e12);
}

real_t EquationOfState::get_rhd_enthalpy(real_t density, real_t pressure, real_t gamma, const std::string& eos_type) {
    if (eos_type == "TM") {
        real_t tau = pressure / density;
        return density * (2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0));
    } else {
        return density + pressure / (gamma - 1.0) + pressure;
    }
}

real_t EquationOfState::get_rhd_sound_speed(real_t density, real_t pressure, real_t gamma, const std::string& eos_type) {
    if (eos_type == "TM") {
        real_t tau = pressure / density;
        real_t enth = 2.5 * tau + 1.5 * std::sqrt(tau*tau + 4.0/9.0);
        return std::sqrt(tau / (3.0 * enth) * (5.0 * enth - 8.0 * tau) / (enth - tau));
    } else {
        real_t enth = density + pressure / (gamma - 1.0) + pressure;
        return std::sqrt(gamma * pressure / enth);
    }
}

real_t EquationOfState::get_rhd_pressure_zero_momentum(real_t D, real_t E, real_t gamma, const std::string& eos_type) {
    if (eos_type == "TM") {
        return (E*E - D*D) / (3.0 * E);
    } else {
        return (E - D) * (gamma - 1.0);
    }
}

real_t EquationOfState::get_rhd_pressure_from_xsi(real_t xsi, real_t rho, real_t gamma, const std::string& eos_type) {
    if (eos_type == "TM") {
        return (2.0 * xsi * (xsi + 2.0 * rho)) / (5.0 * (xsi + rho) + std::sqrt(9.0*xsi*xsi + 18.0*rho*xsi + 25.0*rho*rho));
    } else {
        return (gamma - 1.0) / gamma * xsi;
    }
}

} // namespace ramses
