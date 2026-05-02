#include "ramses/RamsesReader.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Config.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace ramses;

struct CellData {
    double x[3];
    std::vector<double> u;
};

std::vector<CellData> collect_leaf_cells(const AmrGrid& grid) {
    std::vector<CellData> leaf_cells;
    int twotondim = (1 << NDIM);
    real_t dx_base = grid.boxlen; // Simplified
    
    for (int i = 1; i <= grid.ncoarse; ++i) {
        if (grid.son[i-1] == 0) {
            CellData c; 
            // Simplified position
            c.x[0] = (i - 0.5) * (grid.boxlen / grid.ncoarse);
            c.u.resize(grid.nvar);
            for(int iv=1; iv<=grid.nvar; ++iv) c.u[iv-1] = grid.uold(i, iv);
            leaf_cells.push_back(c);
        }
    }

    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int icpu = 1; icpu <= grid.ncpu; ++icpu) {
            int igrid = grid.headl(icpu, il);
            while (igrid > 0) {
                for (int ic = 1; ic <= twotondim; ++ic) {
                    int id = grid.ncoarse + (ic - 1) * grid.ngridmax + igrid;
                    if (grid.son[id - 1] == 0) {
                        CellData c;
                        real_t dx_level = grid.boxlen / (1 << (il - 1));
                        int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                        c.x[0] = grid.xg[0 * grid.ngridmax + (igrid - 1)] + (ix - 0.5) * dx_level;
                        c.x[1] = grid.xg[1 * grid.ngridmax + (igrid - 1)] + (iy - 0.5) * dx_level;
                        c.x[2] = grid.xg[2 * grid.ngridmax + (igrid - 1)] + (iz - 0.5) * dx_level;
                        c.u.resize(grid.nvar);
                        for(int iv=1; iv<=grid.nvar; ++iv) c.u[iv-1] = grid.uold(id, iv);
                        leaf_cells.push_back(c);
                    }
                }
                igrid = grid.next[igrid - 1];
            }
        }
    }
    return leaf_cells;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <ref_amr_file> <loc_amr_file>" << std::endl;
        return 1;
    }

    AmrGrid grid_ref;
    AmrGrid grid_loc;
    RamsesReader reader_ref(argv[1]);
    RamsesReader reader_loc(argv[2]);

    if (!reader_ref.load_amr(grid_ref)) return 1;
    if (!reader_loc.load_amr(grid_loc)) return 1;

    // Load hydro as well
    std::string ref_hydro = argv[1]; ref_hydro.replace(ref_hydro.find("amr"), 3, "hydro");
    std::string loc_hydro = argv[2]; loc_hydro.replace(loc_hydro.find("amr"), 3, "hydro");
    RamsesReader h_ref(ref_hydro); h_ref.load_hydro(grid_ref);
    RamsesReader h_loc(loc_hydro); h_loc.load_hydro(grid_loc);

    auto leaf_ref = collect_leaf_cells(grid_ref);
    auto leaf_loc = collect_leaf_cells(grid_loc);

    std::cout << "Reference leaf cells: " << leaf_ref.size() << std::endl;
    std::cout << "Local leaf cells:     " << leaf_loc.size() << std::endl;

    if (leaf_ref.size() != leaf_loc.size()) {
        std::cout << "Mismatch in leaf cell count!" << std::endl;
        return 1;
    }

    double max_err = 0.0;
    for (size_t i = 0; i < leaf_ref.size(); ++i) {
        for (int iv = 0; iv < grid_ref.nvar; ++iv) {
            double err = std::abs(leaf_ref[i].u[iv] - leaf_loc[i].u[iv]);
            max_err = std::max(max_err, err);
        }
    }

    std::cout << "Maximum absolute error: " << max_err << std::endl;
    if (max_err < 1e-10) {
        std::cout << "Verification PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "Verification FAILED" << std::endl;
        return 1;
    }
}
