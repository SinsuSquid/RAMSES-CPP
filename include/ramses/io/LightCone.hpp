#ifndef RAMSES_LIGHT_CONE_HPP
#define RAMSES_LIGHT_CONE_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Handles light-cone generation for cosmological simulations.
 */
class LightCone {
public:
    LightCone(AmrGrid& grid, Config& config);
    ~LightCone();

    void update(real_t t, real_t dt);
    void output();

private:
    AmrGrid& grid_;
    Config& config_;
    
    real_t z_min_ = 0.0;
    real_t z_max_ = 10.0;
    real_t theta_ = 0.0;
    real_t phi_ = 0.0;
};

} // namespace ramses

#endif // RAMSES_LIGHT_CONE_HPP
