#ifndef RAMSES_TREE_UPDATER_HPP
#define RAMSES_TREE_UPDATER_HPP

#include "AmrGrid.hpp"

namespace ramses {

/**
 * @brief Handles tree refinement and derefinement operations.
 * 
 * Implements logic from amr/refine_utils.f90.
 */
class TreeUpdater {
public:
    TreeUpdater(AmrGrid& grid) : grid_(grid) {}

    /**
     * @brief Refines coarse level cells (Level 1).
     */
    void refine_coarse();

    /**
     * @brief Refines fine level cells (ilevel -> ilevel + 1).
     */
    void refine_fine(int ilevel);

    /**
     * @brief Marks cells for refinement based on gradients.
     */
    void mark_cells(int ilevel);

private:
    AmrGrid& grid_;

    /**
     * @brief Internal helper to create a single grid from a coarse cell.
     */
    void make_grid_coarse(int ind_cell, int ibound, bool boundary_region);

    /**
     * @brief Internal helper to create a single grid from a fine cell.
     */
    void make_grid_fine(int ind_grid_father, int icell_pos, int ilevel, int ibound, bool boundary_region);
};

} // namespace ramses

#endif // RAMSES_TREE_UPDATER_HPP
