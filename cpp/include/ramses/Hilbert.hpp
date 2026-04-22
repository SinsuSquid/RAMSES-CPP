#ifndef RAMSES_HILBERT_HPP
#define RAMSES_HILBERT_HPP

#include "Types.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Hilbert Curve indexing for 1D, 2D, and 3D.
 * 
 * Replicates the state-diagram based implementation from amr/hilbert.f90.
 */
class Hilbert {
public:
    /**
     * @brief Computes 1D Hilbert keys (trivial mapping).
     */
    static void hilbert1d(const std::vector<int>& x, std::vector<qdp_t>& order);

    /**
     * @brief Computes 2D Hilbert keys.
     */
    static void hilbert2d(const std::vector<int>& x, const std::vector<int>& y, 
                         std::vector<qdp_t>& order, int bit_length);

    /**
     * @brief Computes 3D Hilbert keys.
     */
    static void hilbert3d(const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z,
                         std::vector<qdp_t>& order, int bit_length);
};

} // namespace ramses

#endif // RAMSES_HILBERT_HPP
