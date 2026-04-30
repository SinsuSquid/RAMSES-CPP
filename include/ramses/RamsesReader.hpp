#ifndef RAMSES_READER_HPP
#define RAMSES_READER_HPP

#include "AmrGrid.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

namespace ramses {

/**
 * @brief Utility to read RAMSES unformatted Fortran binary files.
 */
class RamsesReader {
public:
    RamsesReader(const std::string& filename) {
        file_.open(filename, std::ios::binary);
    }
    
    ~RamsesReader() {
        if (file_.is_open()) file_.close();
    }

    bool is_open() const { return file_.is_open(); }

    bool load_amr(AmrGrid& grid);
    bool load_hydro(AmrGrid& grid);

private:
    template <typename T>
    T read_single() {
        int32_t size1, size2;
        T val;
        file_.read(reinterpret_cast<char*>(&size1), sizeof(int32_t));
        file_.read(reinterpret_cast<char*>(&val), sizeof(T));
        file_.read(reinterpret_cast<char*>(&size2), sizeof(int32_t));
        return val;
    }

    template <typename T>
    void read_record(std::vector<T>& data) {
        int32_t size1, size2;
        file_.read(reinterpret_cast<char*>(&size1), sizeof(int32_t));
        std::cout << "  [Reader] Reading record of size " << size1 << std::endl;
        if (size1 < 0 || size1 > 1000000000) { std::cerr << "[Reader] Error: Invalid record size " << size1 << std::endl; }
        size_t count = size1 / sizeof(T);
        data.resize(count);
        file_.read(reinterpret_cast<char*>(data.data()), size1);
        file_.read(reinterpret_cast<char*>(&size2), sizeof(int32_t));
    }

    void skip_record() {
        int32_t size1, size2;
        file_.read(reinterpret_cast<char*>(&size1), sizeof(int32_t));
        std::cout << "  [Reader] Skipping record of size " << size1 << std::endl;
        file_.seekg(size1, std::ios::cur);
        file_.read(reinterpret_cast<char*>(&size2), sizeof(int32_t));
    }

    std::ifstream file_;
};

} // namespace ramses

#endif
