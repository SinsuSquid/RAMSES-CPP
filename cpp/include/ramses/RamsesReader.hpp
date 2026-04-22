#ifndef RAMSES_READER_HPP
#define RAMSES_READER_HPP

#include "Types.hpp"
#include "AmrGrid.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

namespace ramses {

class RamsesReader {
public:
    RamsesReader(const std::string& filename) : file_(filename, std::ios::binary) {}
    
    ~RamsesReader() {
        if (file_.is_open()) file_.close();
    }

    bool is_open() const { return file_.is_open(); }

    uint32_t read_tag() {
        uint32_t tag = 0;
        file_.read(reinterpret_cast<char*>(&tag), 4);
        return tag;
    }

    template<typename T>
    size_t read_record(std::vector<T>& data) {
        uint32_t t1 = read_tag();
        if (file_.eof()) return 0;
        
        size_t count = t1 / sizeof(T);
        data.resize(count);
        file_.read(reinterpret_cast<char*>(data.data()), t1);
        
        uint32_t t2 = read_tag();
        if (t1 != t2) {
            std::cerr << "Warning: Fortran record tags do not match! " << t1 << " != " << t2 << std::endl;
        }
        return t1;
    }

    template<typename T>
    T read_single() {
        uint32_t t1 = read_tag();
        T val;
        file_.read(reinterpret_cast<char*>(&val), sizeof(T));
        uint32_t t2 = read_tag();
        return val;
    }

    void skip_record() {
        uint32_t t1 = read_tag();
        file_.seekg(t1, std::ios::cur);
        uint32_t t2 = read_tag();
    }

    bool load_amr(AmrGrid& grid);
    bool load_hydro(AmrGrid& grid);

private:
    std::ifstream file_;
};

} // namespace ramses

#endif // RAMSES_READER_HPP
