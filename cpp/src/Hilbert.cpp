#include "ramses/Hilbert.hpp"
#include <cmath>
#include <iostream>

namespace ramses {

void Hilbert::hilbert1d(const std::vector<int>& x, std::vector<qdp_t>& order) {
    size_t npoint = x.size();
    order.resize(npoint);
    for (size_t ip = 0; ip < npoint; ++ip) {
        order[ip] = static_cast<qdp_t>(x[ip]);
    }
}

void Hilbert::hilbert2d(const std::vector<int>& x, const std::vector<int>& y, 
                       std::vector<qdp_t>& order, int bit_length) {
    size_t npoint = x.size();
    order.assign(npoint, 0.0);

    // Fortran: state_diagram(sdigit,0:1,cstate)
    // sdigit: 0-3, bit: 0-1, cstate: 0-3
    // C++ map: [cstate][bit][sdigit]
    static const int state_diagram[4][2][4] = {
        {{1, 0, 2, 0}, {0, 1, 3, 2}}, // cstate 0
        {{0, 3, 1, 1}, {0, 3, 1, 2}}, // cstate 1
        {{2, 2, 0, 3}, {2, 1, 3, 0}}, // cstate 2
        {{3, 1, 3, 2}, {2, 3, 1, 0}}  // cstate 3
    };

    for (size_t ip = 0; ip < npoint; ++ip) {
        std::vector<bool> i_bit_mask(2 * bit_length);
        
        // Interleave bits
        for (int i = 0; i < bit_length; ++i) {
            i_bit_mask[2 * i + 1] = (x[ip] >> i) & 1;
            i_bit_mask[2 * i]     = (y[ip] >> i) & 1;
        }

        int cstate = 0;
        for (int i = bit_length - 1; i >= 0; --i) {
            int b1 = i_bit_mask[2 * i + 1] ? 1 : 0;
            int b0 = i_bit_mask[2 * i]     ? 1 : 0;
            int sdigit = b1 * 2 + b0;

            int nstate = state_diagram[cstate][0][sdigit];
            int hdigit = state_diagram[cstate][1][sdigit];

            i_bit_mask[2 * i + 1] = (hdigit >> 1) & 1;
            i_bit_mask[2 * i]     = (hdigit >> 0) & 1;
            cstate = nstate;
        }

        order[ip] = 0;
        for (int i = 0; i < 2 * bit_length; ++i) {
            if (i_bit_mask[i]) {
                order[ip] += std::pow(2.0, i);
            }
        }
    }
}

void Hilbert::hilbert3d(const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z,
                       std::vector<qdp_t>& order, int bit_length) {
    size_t npoint = x.size();
    order.assign(npoint, 0.0);

    // Fortran: state_diagram(sdigit,0:1,cstate)
    // sdigit: 0-7, bit: 0-1, cstate: 0-11
    // C++ map: [cstate][bit][sdigit]
    static const int state_diagram[12][2][8] = {
        {{1, 2, 3, 2, 4, 5, 3, 5}, {0, 1, 3, 2, 7, 6, 4, 5}}, // 0
        {{2, 6, 0, 7, 8, 8, 0, 7}, {0, 7, 1, 6, 3, 4, 2, 5}}, // 1
        {{0, 9, 10, 9, 1, 1, 11, 11}, {0, 3, 7, 4, 1, 2, 6, 5}}, // 2
        {{6, 0, 6, 11, 9, 0, 9, 8}, {2, 3, 1, 0, 5, 4, 6, 7}}, // 3
        {{11, 11, 0, 7, 5, 9, 0, 7}, {4, 3, 5, 2, 7, 0, 6, 1}}, // 4
        {{4, 4, 8, 8, 0, 6, 10, 6}, {6, 5, 1, 2, 7, 4, 0, 3}}, // 5
        {{5, 7, 5, 3, 1, 1, 11, 11}, {4, 7, 3, 0, 5, 6, 2, 1}}, // 6
        {{6, 1, 6, 10, 9, 4, 9, 10}, {6, 7, 5, 4, 1, 0, 2, 3}}, // 7
        {{10, 3, 1, 1, 10, 3, 5, 9}, {2, 5, 3, 4, 1, 6, 0, 7}}, // 8
        {{4, 4, 8, 8, 2, 7, 2, 3}, {2, 1, 5, 6, 3, 0, 4, 7}}, // 9
        {{7, 2, 11, 2, 7, 5, 8, 5}, {4, 5, 7, 6, 3, 2, 0, 1}}, // 10
        {{10, 3, 2, 6, 10, 3, 4, 4}, {6, 1, 7, 0, 5, 2, 4, 3}}  // 11
    };

    for (size_t ip = 0; ip < npoint; ++ip) {
        std::vector<bool> i_bit_mask(3 * bit_length);

        for (int i = 0; i < bit_length; ++i) {
            i_bit_mask[3 * i + 2] = (x[ip] >> i) & 1;
            i_bit_mask[3 * i + 1] = (y[ip] >> i) & 1;
            i_bit_mask[3 * i]     = (z[ip] >> i) & 1;
        }

        int cstate = 0;
        for (int i = bit_length - 1; i >= 0; --i) {
            int b2 = i_bit_mask[3 * i + 2] ? 1 : 0;
            int b1 = i_bit_mask[3 * i + 1] ? 1 : 0;
            int b0 = i_bit_mask[3 * i]     ? 1 : 0;
            int sdigit = b2 * 4 + b1 * 2 + b0;

            int nstate = state_diagram[cstate][0][sdigit];
            int hdigit = state_diagram[cstate][1][sdigit];

            i_bit_mask[3 * i + 2] = (hdigit >> 2) & 1;
            i_bit_mask[3 * i + 1] = (hdigit >> 1) & 1;
            i_bit_mask[3 * i]     = (hdigit >> 0) & 1;
            cstate = nstate;
        }

        order[ip] = 0;
        for (int i = 0; i < 3 * bit_length; ++i) {
            if (i_bit_mask[i]) {
                order[ip] += std::pow(2.0, i);
            }
        }
    }
}

} // namespace ramses
