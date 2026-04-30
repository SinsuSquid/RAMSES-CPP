#include "ramses/RamsesReader.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Config.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace ramses;

struct CellData {
    real_t x[3];
    real_t vars[20];
    
    bool operator<(const CellData& other) const {
        if (std::abs(x[0] - other.x[0]) > 1e-12) return x[0] < other.x[0];
        if (std::abs(x[1] - other.x[1]) > 1e-12) return x[1] < other.x[1];
        return x[2] < other.x[2];
    }
};

std::vector<CellData> collect_leaf_cells(const AmrGrid& grid) {
    std::vector<CellData> cells;
    // 1. Coarse cells
    for (int i = 1; i <= grid.ncoarse; ++i) {
        if (grid.son[i] == 0) {
            CellData c; grid.get_cell_center(i, c.x);
            for (int iv = 1; iv <= grid.nvar; ++iv) c.vars[iv - 1] = grid.uold(i, iv);
            cells.push_back(c);
        }
    }
    // 2. Refined cells
    int cells_per_oct = (1 << grid.ndim);
    for (int il = 1; il <= grid.nlevelmax; ++il) {
        for (int ic = 1; ic <= grid.ncpu; ++ic) {
            int ig = grid.headl(ic, il);
            while (ig > 0) {
                for (int icell = 1; icell <= cells_per_oct; ++icell) {
                    int id = grid.ncoarse + (icell - 1) * grid.ngridmax + ig;
                    if (grid.son[id] == 0) {
                        CellData c; grid.get_cell_center(id, c.x);
                        for (int iv = 1; iv <= grid.nvar; ++iv) c.vars[iv - 1] = grid.uold(id, iv);
                        cells.push_back(c);
                    }
                }
                ig = grid.next[ig - 1];
            }
        }
    }
    return cells;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <path_ref_amr> <path_local_amr>" << std::endl;
        return 1;
    }

    std::string ref_amr = argv[1];
    std::string loc_amr = argv[2];
    
    std::string ref_hydro = ref_amr;
    size_t pos = ref_hydro.find("amr_");
    if (pos != std::string::npos) ref_hydro.replace(pos, 4, "hydro_");

    std::string loc_hydro = loc_amr;
    pos = loc_hydro.find("amr_");
    if (pos != std::string::npos) loc_hydro.replace(pos, 4, "hydro_");

    Config config;
    AmrGrid grid_ref(config);
    AmrGrid grid_loc(config);

    std::cout << "Loading Reference..." << std::endl;
    {
        RamsesReader reader(ref_amr);
        if (!reader.load_amr(grid_ref)) return 1;
        RamsesReader hreader(ref_hydro);
        if (!hreader.load_hydro(grid_ref)) return 1;
    }

    std::cout << "Loading Local..." << std::endl;
    {
        RamsesReader reader(loc_amr);
        if (!reader.load_amr(grid_loc)) return 1;
        RamsesReader hreader(loc_hydro);
        if (!hreader.load_hydro(grid_loc)) return 1;
    }

    std::vector<CellData> cells_ref = collect_leaf_cells(grid_ref);
    std::vector<CellData> cells_loc = collect_leaf_cells(grid_loc);

    if (cells_ref.size() != cells_loc.size()) {
        std::cerr << "Leaf cell count mismatch: " << cells_ref.size() << " vs " << cells_loc.size() << std::endl;
        return 1;
    }

    std::cout << "Sorting " << cells_ref.size() << " cells..." << std::endl;
    std::sort(cells_ref.begin(), cells_ref.end());
    std::sort(cells_loc.begin(), cells_loc.end());

    std::cout << "Comparing primitive variables..." << std::endl;
    bool all_ok = true;
    int n_comp = 1 + grid_ref.ndim + 1;
    for (int iv = 0; iv < n_comp; ++iv) {
        real_t max_diff = 0;
        real_t avg_diff = 0;
        for (size_t i = 0; i < cells_ref.size(); ++i) {
            real_t diff = std::abs(cells_ref[i].vars[iv] - cells_loc[i].vars[iv]);
            max_diff = std::max(max_diff, diff);
            avg_diff += diff;
        }
        avg_diff /= cells_ref.size();
        std::string var_name = (iv == 0) ? "Density" : (iv <= grid_ref.ndim ? "Velocity" : "Pressure");
        std::cout << "Variable " << iv+1 << " (" << var_name << "): Max Diff = " << max_diff << ", Avg Diff = " << avg_diff << std::endl;
        if (max_diff > 1e-8) all_ok = false;
    }

    if (all_ok) std::cout << "\nSUCCESS: Snapshots match within tolerance." << std::endl;
    else std::cout << "\nFAILURE: Discrepancies detected." << std::endl;

    return all_ok ? 0 : 1;
}
