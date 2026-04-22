#ifndef RAMSES_CONSTANTS_HPP
#define RAMSES_CONSTANTS_HPP

#include "Types.hpp"

namespace ramses {
namespace constants {

#if NDIM == 1
    constexpr int twotondim = 2;
    constexpr int threetondim = 3;
#elif NDIM == 2
    constexpr int twotondim = 4;
    constexpr int threetondim = 9;
#elif NDIM == 3
    constexpr int twotondim = 8;
    constexpr int threetondim = 27;
#endif
    constexpr int twondim = 2 * NDIM;
    constexpr int MAXBOUND = 100;

    // Lookup tables from amr_constants.f90 (iii, jjj)
    // Layout: [ndim][2][twotondim]
    extern const int iii[3][2][8];
    extern const int jjj[3][2][8];

    // Lookup tables for 3x3x3 gathering (lll, mmm)
    // Layout: [twotondim][threetondim]
    extern const int lll[8][27];
    extern const int mmm[8][27];
    
} // namespace constants
} // namespace ramses

#endif // RAMSES_CONSTANTS_HPP
