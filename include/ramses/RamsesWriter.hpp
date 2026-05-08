#ifndef RAMSES_WRITER_HPP
#define RAMSES_WRITER_HPP

#include "AmrGrid.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>

namespace ramses {

/**
 * @brief Metadata for a simulation snapshot.
 */
struct SnapshotInfo {
    real_t t;
    int nstep;
    int nstep_coarse;
    int noutput;
    int iout;
    std::vector<double> tout;
    real_t gamma;
    int nener = 0;
    
    // Legacy fields
    real_t einit = 0.0, mass_tot_0 = 0.0, rho_tot = 0.0;
    real_t omega_m = 0.0, omega_l = 0.0, omega_k = 0.0, omega_b = 0.0;
    real_t h0 = 0.0, aexp_ini = 1.0, boxlen_ini = 1.0;
    real_t aexp = 1.0, hexp = 0.0, aexp_old = 1.0;
    real_t epot_tot_int = 0.0, epot_tot_old = 0.0;
    real_t mass_sph = 0.0;
};

class RamsesWriter {
public:
    RamsesWriter(const std::string& filename);
    bool is_open() const;
    void write_amr(const AmrGrid& grid, const SnapshotInfo& info);
    void write_hydro(const AmrGrid& grid, const SnapshotInfo& info);
    void write_grav(const AmrGrid& grid, const SnapshotInfo& info);
    void write_header(const AmrGrid& grid, const SnapshotInfo& info);
    void write_header_file(const AmrGrid& grid, const SnapshotInfo& info);
    void write_hydro_descriptor(const AmrGrid& grid, const SnapshotInfo& info);
    void write_rt(const AmrGrid& grid, const SnapshotInfo& info, int nGroups, real_t rt_c);
    void write_rt_descriptor(const AmrGrid& grid, const SnapshotInfo& info, int nGroups);
    void write_particles(const AmrGrid& grid, const SnapshotInfo& info);
    void write_particles_descriptor(const AmrGrid& grid, const SnapshotInfo& info);

private:
    template <typename T>
    void write_record_internal(std::ofstream& file, const T* data, size_t count) {
        int32_t size = count * sizeof(T);
        file.write((char*)&size, 4);
        if (count > 0) file.write((char*)data, size);
        file.write((char*)&size, 4);
    }
    
    template <typename T>
    void write_record(const T* data, size_t count) {
        std::ofstream file(filename_, std::ios::binary | std::ios::app);
        write_record_internal(file, data, count);
    }
    
    std::string filename_;
};

} // namespace ramses

#endif
