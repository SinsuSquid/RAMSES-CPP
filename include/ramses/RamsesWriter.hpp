#ifndef RAMSES_WRITER_HPP
#define RAMSES_WRITER_HPP

#include "AmrGrid.hpp"
#include <fstream>
#include <string>
#include <cstdint>

namespace ramses {

/**
 * @brief Metadata for a simulation snapshot.
 */
struct SnapshotInfo {
    real_t t;
    int nstep;
    int noutput;
    int iout;
    std::vector<double> tout;
    real_t gamma;
    int nener = 0;
};

class RamsesWriter {
public:
    RamsesWriter(const std::string& filename);
    bool is_open() const;
    void write_amr(const AmrGrid& grid, const SnapshotInfo& info);
    void write_hydro(const AmrGrid& grid, const SnapshotInfo& info);
    void write_grav(const AmrGrid& grid, const SnapshotInfo& info);

private:
    template <typename T>
    void write_record(const T* data, size_t count);
    std::ofstream file_;
};

} // namespace ramses

#endif
