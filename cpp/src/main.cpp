#include <iostream>
#include <vector>
#include "ramses/AmrGrid.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/TreeUpdater.hpp"

int main() {
    std::cout << "--- Testing Tree Refinement ---" << std::endl;
    
    namespace p = ramses::params;
    p::nx = 2; p::ny = 2; p::nz = 2;
    int ngridmax = 100;
    int nvar = 1;
    int ncpu = 1;
    int nlevelmax = 10;
    
    ramses::AmrGrid grid;
    grid.allocate(p::nx, p::ny, p::nz, ngridmax, nvar, ncpu, nlevelmax);
    
    // Mark cell (1,1,1) for refinement. 
    grid.flag1[1] = 1;
    grid.flag2[1] = 1;
    grid.cpu_map[1] = 1; // 1-based icpu
    
    ramses::TreeUpdater updater(grid);
    updater.refine_coarse();
    
    std::cout << "Refinement complete." << std::endl;
    
    if (grid.son[1] > 0) {
        int igrid = grid.son[1];
        std::cout << "Cell 1 refined. New Grid Index: " << igrid << std::endl;
        std::cout << "Grid " << igrid << " Center: (" 
                  << grid.get_xg(igrid, 1) << ", "
                  << grid.get_xg(igrid, 2) << ", "
                  << grid.get_xg(igrid, 3) << ")" << std::endl;
        
        std::cout << "Grid " << igrid << " Neighbors: ";
        for (int i = 1; i <= 6; ++i) {
            std::cout << grid.get_nbor(igrid, i) << " ";
        }
        std::cout << std::endl;
        
        // Neighbor 2 should be cell (2,1,1) -> ix=1, iy=0, iz=0 -> ind = 2
        if (grid.get_nbor(igrid, 2) == 2) {
            std::cout << "Right neighbor correct (Cell 2)." << std::endl;
        }
    } else {
        std::cout << "Refinement FAILED." << std::endl;
        return 1;
    }

    return 0;
}
