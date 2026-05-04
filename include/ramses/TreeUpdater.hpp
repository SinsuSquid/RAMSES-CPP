#ifndef RAMSES_TREE_UPDATER_HPP
#define RAMSES_TREE_UPDATER_HPP

#include "AmrGrid.hpp"
#include "Config.hpp"
#include <functional>

namespace ramses {

/**
 * @brief Handles tree refinement and derefinement operations.
 */
class TreeUpdater {
public:
    TreeUpdater(AmrGrid& grid, Config& config);

    void make_grid_fine(int ilevel);
    void remove_grid_fine(int ilevel);
    void restrict_fine(int ilevel);
    void flag_fine(int ilevel, real_t err_grad_d, real_t err_grad_p, real_t err_grad_v, real_t err_grad_b2 = -1.0);

    using InterpolHook = std::function<void(const real_t u1[7][20], real_t u2[8][20])>;
    void set_interpol_hook(InterpolHook hook) { interpol_hook_ = hook; }

private:
    AmrGrid& grid_;
    Config& config_;
    InterpolHook interpol_hook_;
};

} // namespace ramses

#endif // RAMSES_TREE_UPDATER_HPP
