#ifndef RAMSES_WRITER_HPP
#define RAMSES_WRITER_HPP

#include "Types.hpp"
#include "AmrGrid.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace ramses {

/**
 * @brief Writes RAMSES unformatted Fortran binary files.
 */
class RamsesWriter {
public:
    RamsesWriter(const std::string& filename) : file_(filename, std::ios::binary) {}
    ~RamsesWriter() { if (file_.is_open()) file_.close(); }

    bool is_open() const { return file_.is_open(); }

    void write_tag(uint32_t tag) {
        file_.write(reinterpret_cast<const char*>(&tag), 4);
    }

    template<typename T>
    void write_record(const std::vector<T>& data) {
        uint32_t tag = static_cast<uint32_t>(data.size() * sizeof(T));
        write_tag(tag);
        file_.write(reinterpret_cast<const char*>(data.data()), tag);
        write_tag(tag);
    }

    template<typename T>
    void write_single(T val) {
        uint32_t tag = sizeof(T);
        write_tag(tag);
        file_.write(reinterpret_cast<const char*>(&val), tag);
        write_tag(tag);
    }

    /**
     * @brief Dumps the full AMR grid to a binary file.
     */
    void write_amr(const AmrGrid& grid);

    /**
     * @brief Dumps the hydro variables to a binary file.
     */
    void write_hydro(const AmrGrid& grid, int ilevel);

private:
    std::ofstream file_;
};

} // namespace ramses

#endif // RAMSES_WRITER_HPP
