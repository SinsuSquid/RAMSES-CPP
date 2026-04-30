#ifndef RAMSES_AMR_GRID_HPP
#define RAMSES_AMR_GRID_HPP

#include "Types.hpp"
#include "Field.hpp"
#include "Constants.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Manages the AMR tree structure and associated physical fields.
 * 
 * Replicates the data structures in amr_commons.f90 and hydro_commons.f90.
 */
class AmrGrid {
public:
    AmrGrid() = default;

    /**
     * @brief Allocates memory for the grid based on dimensions and limits.
     */
    void allocate(int nx, int ny, int nz, int ngridmax, int nvar, int ncpu, int nlevelmax);

    // Tree Structure (Oct-based)
    std::vector<real_t> xg;       // Grid positions [ngridmax * NDIM]
    std::vector<int> father;      // Father cell index [ngridmax]
    std::vector<int> nbor;        // Neighboring father cells [ngridmax * 2*NDIM]
    std::vector<int> next;        // Next grid in linked list [ngridmax]
    std::vector<int> prev;        // Previous grid in linked list [ngridmax]

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
    Field<real_t> f;              // Gravitational acceleration [ncell, NDIM]
    std::vector<real_t> rho;      // Density (gravity source) [ncell]

    // Linked list pointers (1-based level indexing)
    Field<int> headl;             // Head grid in level list [ncpu, nlevelmax]
    Field<int> taill;             // Tail grid in level list [ncpu, nlevelmax]
    Field<int> numbl;             // Number of grids in level list [ncpu, nlevelmax]

    Field<int> headb;             // Head grid in boundary list [MAXBOUND, nlevelmax]
    Field<int> tailb;             // Tail grid in boundary list [MAXBOUND, nlevelmax]
    Field<int> numbb;             // Number of grids in boundary list [MAXBOUND, nlevelmax]

    // Linked list pointers for free memory

    int headf, tailf, numbf;      // Free memory list pointers
    
    // Grid info
    int ncoarse;
    int ncell;
    int ngridmax;
    int nvar;
    int ncpu;
    int nlevelmax;
    int ndim;
    real_t rho_tot = 0.0;

    /**
     * @brief Returns total number of grids at a specific level across all CPUs.
     */
    int count_grids_at_level(int ilevel) const {
        int total = 0;
        for (int i = 1; i <= ncpu; ++i) {
            total += numbl(i, ilevel);
        }
        return total;
    }

    // 1-based utility accessors
    inline real_t& get_xg(int igrid, int idim) { return xg[(idim - 1) * ngridmax + (igrid - 1)]; }
    inline int& get_nbor(int igrid, int iface) { return nbor[(iface - 1) * ngridmax + (igrid - 1)]; }

    /**
     * @brief Find neighboring grids of an oct.
     * @param igrid 1-based oct index.
     * @param igridn Output 1-based indices of neighboring grids (0 if not refined). igridn[0] is center.
     */
    void get_nbor_grids(int igrid, int igridn[7]) const;

    /**
     * @brief Find neighboring cells of a cell within an oct.
     * @param igridn 1-based indices of neighboring grids (from get_nbor_grids).
     * @param icell_pos 1-based position within oct (1-8).
     * @param icelln Output 1-based indices of neighboring cells.
     */
    void get_nbor_cells(const int igridn[7], int icell_pos, int icelln[6]) const;

    /**
     * @brief Find the 3x3x3 cube of father cells around an oct.
     * @param igrid 1-based oct index.
     * @param nbors_father Output array of 27 cell indices.
     */
    void get_3x3x3_father(int igrid, int nbors_father[27]) const;

    /**
     * @brief Computes the center coordinates of a cell.
     * @param icell 1-based cell index.
     * @param x Output center coordinates.
     */
    void get_cell_center(int icell, real_t x[3]) const;

private:
    // Helper to calculate ncell = ncoarse + twotondim * ngridmax
    static int calculate_ncell(int nx, int ny, int nz, int ngridmax) {
        return (nx * ny * nz) + constants::twotondim * ngridmax;
    }
};

} // namespace ramses

#endif // RAMSES_AMR_GRID_HPP
