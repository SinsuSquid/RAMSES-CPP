#ifndef RAMSES_TYPES_HPP
#define RAMSES_TYPES_HPP

#include <cstdint>

namespace ramses {

// Precision types mapping to Fortran sp, dp, qdp
#if RAMSES_PRECISION == 4
    using real_t = float;
#else
    using real_t = double;
#endif

using sp_t = float;
using dp_t = double;
using qdp_t = double; // Or __float128 if QUADHILBERT is needed

// Integer types mapping to i4b, i8b
using i4b_t = int32_t;

#ifdef LONGINT
    using i8b_t = int64_t;
#else
    using i8b_t = int32_t;
#endif

} // namespace ramses

#endif // RAMSES_TYPES_HPP
