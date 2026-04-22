#ifndef RAMSES_FIELD_HPP
#define RAMSES_FIELD_HPP

#include "Types.hpp"
#include <vector>
#include <stdexcept>

namespace ramses {

/**
 * @brief A 2D array wrapper designed to mimic Fortran's (1:N, 1:M) indexing.
 * 
 * Uses a flat std::vector for memory management.
 */
template<typename T>
class Field {
public:
    Field() : dim1_(0), dim2_(0) {}
    Field(size_t n, size_t m) : dim1_(n), dim2_(m), data_(n * m, T(0)) {}

    void allocate(size_t n, size_t m) {
        dim1_ = n;
        dim2_ = m;
        data_.assign(n * m, T(0));
    }

    // 1-based access: field(i, j)
    inline T& operator()(size_t i, size_t j) {
        return data_[(j - 1) * dim1_ + (i - 1)];
    }

    inline const T& operator()(size_t i, size_t j) const {
        return data_[(j - 1) * dim1_ + (i - 1)];
    }

    // Accessors for dimensions
    size_t dim1() const { return dim1_; }
    size_t dim2() const { return dim2_; }

    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

private:
    size_t dim1_, dim2_;
    std::vector<T> data_;
};

} // namespace ramses

#endif // RAMSES_FIELD_HPP
