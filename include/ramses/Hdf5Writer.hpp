#ifndef RAMSES_HDF5_WRITER_HPP
#define RAMSES_HDF5_WRITER_HPP

#include "AmrGrid.hpp"
#include <string>

namespace ramses {

struct SnapshotInfo; // Forward declaration

/**
 * @brief Writer for HDF5 formatted snapshots.
 */
class Hdf5Writer {
public:
    Hdf5Writer(const std::string& filename);
    ~Hdf5Writer();

    void write_snapshot(const AmrGrid& grid, const SnapshotInfo& info);

private:
    std::string filename_;
};

} // namespace ramses

#endif // RAMSES_HDF5_WRITER_HPP
