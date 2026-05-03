#include "ramses/AmrGrid.hpp"
#include "ramses/RamsesReader.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace ramses;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <ref_amr> <local_amr>" << std::endl;
        return 1;
    }

    std::string ref_amr = argv[1];
    std::string loc_amr = argv[2];

    auto get_hydro_path = [](const std::string& amr_path) {
        std::string hydro_path = amr_path;
        size_t pos = hydro_path.find("amr");
        if (pos != std::string::npos) {
            hydro_path.replace(pos, 3, "hydro");
        }
        return hydro_path;
    };

    std::string ref_hydro = get_hydro_path(ref_amr);
    std::string loc_hydro = get_hydro_path(loc_amr);

    AmrGrid grid_ref, grid_loc;
    grid_ref.nvar = 5; grid_loc.nvar = 5; // Default for check

    RamsesReader reader_ref(ref_amr);
    if (!reader_ref.load_amr(grid_ref)) {
        std::cerr << "Failed to load ref AMR: " << ref_amr << std::endl;
        return 1;
    }
    RamsesReader reader_ref_h(ref_hydro);
    reader_ref_h.load_hydro(grid_ref);

    RamsesReader reader_loc(loc_amr);
    if (!reader_loc.load_amr(grid_loc)) {
        std::cerr << "Failed to load local AMR: " << loc_amr << std::endl;
        return 1;
    }
    RamsesReader reader_loc_h(loc_hydro);
    reader_loc_h.load_hydro(grid_loc);

    std::cout << "Comparing " << ref_amr << " and " << loc_amr << std::endl;
    
    // Simple comparison of cell counts and some values
    if (grid_ref.ncoarse != grid_loc.ncoarse) {
        std::cout << "DIFFERENCE: ncoarse ref=" << grid_ref.ncoarse << " loc=" << grid_loc.ncoarse << std::endl;
    }
    
    int diff_count = 0;
    for (int i = 0; i < std::min(grid_ref.ncell, grid_loc.ncell); ++i) {
        for (int iv = 1; iv <= std::min(grid_ref.nvar, grid_loc.nvar); ++iv) {
            if (std::abs(grid_ref.uold(i + 1, iv) - grid_loc.uold(i + 1, iv)) > 1e-10) {
                if (diff_count < 10) {
                    std::cout << "Value diff at cell " << i+1 << " var " << iv << ": ref=" << grid_ref.uold(i+1, iv) << " loc=" << grid_loc.uold(i+1, iv) << std::endl;
                }
                diff_count++;
            }
        }
    }

    if (diff_count == 0 && grid_ref.ncell == grid_loc.ncell) {
        std::cout << "SUCCESS: Snapshots match perfectly!" << std::endl;
        return 0;
    } else {
        std::cout << "FAILURE: Found " << diff_count << " differences." << std::endl;
        return 1;
    }
}
