#include "ramses/LightCone.hpp"
#include <iostream>

namespace ramses {

LightCone::LightCone(AmrGrid& grid, Config& config) : grid_(grid), config_(config) {
    z_min_ = config_.get_double("cosmo_params", "z_min", 0.0);
    z_max_ = config_.get_double("cosmo_params", "z_max", 10.0);
}

LightCone::~LightCone() {}

void LightCone::update(real_t t, real_t dt) {
    // Logic to identify cells/particles crossing the light-cone shell
}

void LightCone::output() {
    // Logic to write light-cone snapshots
}

} // namespace ramses
