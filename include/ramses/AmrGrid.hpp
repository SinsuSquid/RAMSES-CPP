#ifndef RAMSES_AMR_GRID_HPP
#define RAMSES_AMR_GRID_HPP

#include "Types.hpp"
#include "Field.hpp"
#include "Constants.hpp"
#include "Config.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Manages the AMR tree structure and associated physical fields.
 * 
 * Replicates the data structures in amr_commons.f90 and hydro_commons.f90.
 */
class AmrGrid {
public:
    AmrGrid(Config& config) : config_(config) {}

    Config& config() { return config_; }

    /**
     * @brief Allocates memory for the grid based on dimensions and limits.
     */
    void allocate(int nx, int ny, int nz, int ngridmax, int nvar, int ncpu, int nlevelmax);

    // Tree Structure (Oct-based)
    std::vector<real_t> xg;       // Grid center [ngridmax * NDIM]
    std::vector<int> father;      // Father cell index [ngridmax]
    std::vector<int> nbor;        // Neighbor cell indices [ngridmax * 2 * NDIM]
    std::vector<int> next;        // Next grid in level list [ngridmax]
    std::vector<int> prev;        // Previous grid in level list [ngridmax]

    // Cell-based arrays
    std::vector<int> son;         // Son grid index [ncell]
    std::vector<int> flag1;       // Refinement flag [ncell]
    std::vector<int> flag2;       // Expansion flag [ncell]
    std::vector<int> cpu_map;     // Domain decomposition [ncell]
    std::vector<qdp_t> hilbert_keys; // Hilbert keys for load balancing [ncell]
    
    // Physical Fields
    Field<real_t> uold;           // State vector [ncell, nvar]
    Field<real_t> unew;           // Updated state vector [ncell, nvar]
    std::vector<real_t> divu;     // Velocity divergence [ncell]
    std::vector<real_t> phi;      // Gravitational potential [ncell]
    Field<real_t> f;              // Gravitational force [ncell, NDIM]
    std::vector<real_t> rho;      // Density [ncell]

    // Linked List Pointers
    Field<int> headl;             // Head grid in level list [ncpu, nlevelmax]
    Field<int> taill;             // Tail grid in level list [ncpu, nlevelmax]
    Field<int> numbl;             // Number of grids in level list [ncpu, nlevelmax]

    Field<int> headb;             // Head grid in boundary list [MAXBOUND, nlevelmax]
    Field<int> tailb;             // Tail grid in boundary list [MAXBOUND, nlevelmax]
    Field<int> numbb;             // Number of grids in boundary list [MAXBOUND, nlevelmax]

    // Free list management
    int headf, tailf, numbf;
    
    // Metadata
    int ncoarse = 0, ngridmax = 0, nvar = 0, ncpu = 0, nlevelmax = 0, ndim = 0, ncell = 0;

    // Helper to find neighboring grids and cells
    void get_nbor_grids(int igrid, int igridn[7]) const;
    void get_nbor_cells(const int igridn[7], int icell_pos, int icelln[6], int igrid) const;
    void get_27_cell_neighbors(int icell, int nbors[27]) const;

    void restrict_coarse();
    void restrict_fine(int ilevel);
    
    void setup_root_periodicity();

    /**
     * @brief Computes the number of grids currently allocated at a specific level.
     */
    int count_grids_at_level(int ilevel) const {
        if (ilevel <= 0 || ilevel > nlevelmax) return 0;
        int count = 0;
        for (int i = 1; i <= ncpu; ++i) count += numbl(i, ilevel);
        return count;
    }
    
    /**
     * @brief Returns the head of the linked list for a given CPU and level.
     */
    int get_headl(int icpu, int ilevel) const {
        if (ilevel <= 0 || ilevel > nlevelmax) return 0;
        return headl(icpu, ilevel);
    }

    /**
     * @brief Computes the center coordinates of a given cell.
     * @param icell 1-based index of the cell.
     * @param x Output center coordinates.
     */
    void get_cell_center(int icell, real_t x[3]) const;

    /**
     * @brief Find a cell by its normalized coordinates [0, 1].
     */
    int find_cell_by_coords(const real_t x[3], int ilevel_max = -1) const;

    inline real_t& get_xg(int igrid, int idim) { return xg[(idim - 1) * ngridmax + (igrid - 1)]; }
    inline const real_t& get_xg(int igrid, int idim) const { return xg[(idim - 1) * ngridmax + (igrid - 1)]; }
    inline int& get_nbor(int igrid, int inbor) { return nbor[(inbor - 1) * ngridmax + (igrid - 1)]; }
    inline const int& get_nbor(int igrid, int inbor) const { return nbor[(inbor - 1) * ngridmax + (igrid - 1)]; }

    real_t gamma = 1.4;
    real_t err_grad_d = 0.05;
    real_t err_grad_u = 0.05;
    real_t err_grad_p = 0.05;

private:
    // Helper to calculate ncell = ncoarse + twotondim * ngridmax
    static int calculate_ncell(int nx, int ny, int nz, int ngridmax) {
        return (nx * ny * nz) + constants::twotondim * ngridmax;
    }

    Config& config_;
};

} // namespace ramses

#endif // RAMSES_AMR_GRID_HPP
