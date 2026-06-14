#include "ramses/particles/ParticleSystem.hpp"
#include "ramses/core/Constants.hpp"

namespace ramses {

void ParticleSystem::allocate(int npartmax_val, int ngridmax_val) {
    npartmax = npartmax_val;
    ngridmax = ngridmax_val;

    // Allocate particle data
    xp.assign(static_cast<size_t>(npartmax) * NDIM, 0.0);
    vp.assign(static_cast<size_t>(npartmax) * NDIM, 0.0);
    mp.assign(npartmax, 0.0);
    levelp.assign(npartmax, 0);
    idp.assign(npartmax, 0);

    // Linked list pointers
    nextp.assign(npartmax, 0);
    prevp.assign(npartmax, 0);

    // Grid pointers (sized for coarse cells + fine grids)
    // Actually in RAMSES it is dimensioned by (1:npartmax) for headp/tailp? 
    // No, headp is dimensioned by (1:ncell)
    // We'll use ncell-like size here.
    size_t ncell_ptr = 1 + ngridmax * constants::twotondim; // Simplified
    headp.assign(ncell_ptr, 0);
    tailp.assign(ncell_ptr, 0);
    numbp.assign(ncell_ptr, 0);

    // Initialize free list
    for (int i = 1; i < npartmax; ++i) {
        nextp[i - 1] = i + 1;
    }
    for (int i = 2; i <= npartmax; ++i) {
        prevp[i - 1] = i - 1;
    }
    
    headp_free = 1;
    tailp_free = npartmax;
    numbp_free = npartmax;
    
    if (npartmax > 0) {
        prevp[headp_free - 1] = 0;
        nextp[tailp_free - 1] = 0;
    }
}

} // namespace ramses
