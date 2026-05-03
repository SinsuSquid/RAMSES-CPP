#ifndef RAMSES_READER_HPP
#define RAMSES_READER_HPP

#include "AmrGrid.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace ramses {

class RamsesReader {
public:
    RamsesReader(const std::string& filename);
    bool load_amr(AmrGrid& grid);
    bool load_hydro(AmrGrid& grid);

private:
    template <typename T> void read_record(std::ifstream& f, T& data);
    template <typename T> void read_record(std::ifstream& f, std::vector<T>& data);
    void read_record(std::ifstream& f, void* data, size_t size);
    void read_record_raw(std::ifstream& f, void* data, size_t size);

    std::string filename_;
};

} // namespace ramses
#endif
