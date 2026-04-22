#include "ramses/Constants.hpp"

namespace ramses {
namespace constants {

// iii[ndim][2][twotondim]
// Replicated from amr_constants.f90
const int iii[3][2][8] = {
    // idim = 0 (x)
    {
        {1,0,1,0,1,0,1,0}, // inbor = 1
        {0,2,0,2,0,2,0,2}  // inbor = 2
    },
    // idim = 1 (y)
    {
        {3,3,0,0,3,3,0,0}, // inbor = 1
        {0,0,4,4,0,0,4,4}  // inbor = 2
    },
    // idim = 2 (z)
    {
        {5,5,5,5,0,0,0,0}, // inbor = 1
        {0,0,0,0,6,6,6,6}  // inbor = 2
    }
};

// jjj[ndim][2][twotondim]
const int jjj[3][2][8] = {
    // idim = 0 (x)
    {
        {2,1,4,3,6,5,8,7},
        {2,1,4,3,6,5,8,7}
    },
    // idim = 1 (y)
    {
        {3,4,1,2,7,8,5,6},
        {3,4,1,2,7,8,5,6}
    },
    // idim = 2 (z)
    {
        {5,6,7,8,1,2,3,4},
        {5,6,7,8,1,2,3,4}
    }
};

} // namespace constants
} // namespace ramses
