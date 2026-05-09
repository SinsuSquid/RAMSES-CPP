#ifndef RAMSES_CONSTANTS_HPP
#define RAMSES_CONSTANTS_HPP

#include "Types.hpp"

namespace ramses {
namespace constants {

#if NDIM == 1
    constexpr int twotondim = 2;
#elif NDIM == 2
    constexpr int twotondim = 4;
#else
    constexpr int twotondim = 8;
#endif

    constexpr int twondim = 2 * NDIM;
    constexpr int threetondim = (NDIM == 1) ? 3 : ((NDIM == 2) ? 9 : 27);

    // Neighbor and Stencil lookup tables
    extern const int iii[3][2][8];
    extern const int jjj[3][2][8];
    extern const int lll[8][27];
    extern const int mmm[8][27];

    // Coordinate offsets
    extern const real_t xcent[8][3];

} // namespace constants
} // namespace ramses

#endif // RAMSES_CONSTANTS_HPP
