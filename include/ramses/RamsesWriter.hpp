#ifndef RAMSES_WRITER_HPP
#define RAMSES_WRITER_HPP

#include "AmrGrid.hpp"
#include <fstream>
#include <string>
#include <cstdint>

namespace ramses {

class RamsesWriter {
public:
    RamsesWriter(const std::string& filename);
    bool is_open() const;
    void write_amr(const AmrGrid& grid);
    void write_hydro(const AmrGrid& grid, int nlevelmax);

private:
    template <typename T>
    void write_record(const T* data, size_t count);
    std::ofstream file_;
};

} // namespace ramses

#endif
