#include "ramses/Initializer.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>

namespace ramses {

Initializer::~Initializer() {}

void Initializer::apply_all() {
    namespace p = ramses::params;
    std::string filetype = config_.get("init_params", "filetype", "region");
    if (filetype == "grafic") {
        load_grafic();
        return;
    }

    std::string kind = config_.get("init_params", "condinit_kind", "region");
    if (kind == "orzag_tang") {
        real_t pi = std::acos(-1.0);
        real_t B0 = 1.0 / std::sqrt(4.0 * pi);
        real_t d0 = 25.0 / (36.0 * pi);
        real_t p0 = 5.0 / (12.0 * pi);
        real_t gamma = grid_.gamma;

        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int head = grid_.get_headl(1, il);
            int igrid = head;
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << (il - 1)));
            while (igrid > 0) {
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid - 1;
                    int ixyz[3]; int idx = ic - 1; ixyz[2] = idx / 4; idx %= 4; ixyz[1] = idx / 2; ixyz[0] = idx % 2;
                    real_t xc = grid_.xg[0 * grid_.ngridmax + igrid - 1] + (ixyz[0] - 0.5) * dx;
                    real_t yc = (NDIM > 1) ? grid_.xg[1 * grid_.ngridmax + igrid - 1] + (ixyz[1] - 0.5) * dx : 0.5;
                    
                    real_t vx = -std::sin(2.0 * pi * yc);
                    real_t vy = (NDIM > 1) ? std::sin(2.0 * pi * xc) : 0.0;
                    
                    grid_.uold(idc + 1, 1) = d0;
                    grid_.uold(idc + 1, 2) = d0 * vx;
                    grid_.uold(idc + 1, 3) = d0 * vy;
                    grid_.uold(idc + 1, 4) = 0.0;

                    // Magnetic field (face-centered)
                    real_t xl = xc - 0.5 * dx, xr = xc + 0.5 * dx;
                    real_t yl = yc - 0.5 * dx, yr = yc + 0.5 * dx;
                    auto A = [&](real_t x, real_t y) { return B0 * (std::cos(4.0 * pi * x) / (4.0 * pi) + std::cos(2.0 * pi * y) / (2.0 * pi)); };
                    
                    grid_.uold(idc + 1, 6) = (A(xc, yr) - A(xc, yl)) / dx; // Bx_left
                    grid_.uold(idc + 1, grid_.nvar - 2) = (A(xc, yr) - A(xc, yl)) / dx; // Bx_right (simplified)
                    grid_.uold(idc + 1, 7) = (A(xl, yc) - A(xr, yc)) / dx; // By_left
                    grid_.uold(idc + 1, grid_.nvar - 1) = (A(xl, yc) - A(xr, yc)) / dx; // By_right

                    real_t em = 0.5 * (std::pow(grid_.uold(idc+1, 6), 2) + std::pow(grid_.uold(idc+1, 7), 2));
                    grid_.uold(idc + 1, 5) = p0 / (gamma - 1.0) + 0.5 * d0 * (vx * vx + vy * vy) + em;
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    } else if (kind == "ana_disk_potential") {
        real_t dens0 = 0.66;
        real_t height0 = 150.0;
        real_t temp0 = 8000.0;
        real_t gamma = grid_.gamma;
        real_t mu_gas = config_.get_double("cooling_params", "mu_gas", 1.4);
        real_t scale_v = config_.get_double("units_params", "units_length", 1.0) / config_.get_double("units_params", "units_time", 1.0);
        real_t kB = 1.3806e-16, mH = 1.6605e-24;

        auto apply_ana = [&](int idc, real_t x, real_t y, real_t z) {
            real_t rz = (NDIM == 1 ? x : (NDIM == 2 ? y : z)) - 0.5 * params::boxlen;
            real_t d = dens0 * std::exp(-(rz * rz) / (2.0 * height0 * height0));
            grid_.uold(idc, 1) = d;
            grid_.uold(idc, 2) = 0.0; grid_.uold(idc, 3) = 0.0; grid_.uold(idc, 4) = 0.0;
            real_t p = d * (kB * temp0 / (mu_gas * mH)) / (scale_v * scale_v);
            grid_.uold(idc, 5) = p / (gamma - 1.0);
        };

        for (int i = 1; i <= grid_.ncoarse; ++i) {
            int idx = i - 1; int iz = idx / (grid_.nx * grid_.ny); idx %= (grid_.nx * grid_.ny); int iy = idx / grid_.nx; int ix = idx % grid_.nx;
            real_t dx_coarse = grid_.boxlen / std::max({grid_.nx, grid_.ny, grid_.nz});
            apply_ana(i, (ix + 0.5) * dx_coarse, (iy + 0.5) * dx_coarse, (iz + 0.5) * dx_coarse);
        }
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            int igrid = grid_.get_headl(1, il);
            real_t dx = params::boxlen / (real_t)(params::nx * (1 << (il - 1)));
            while (igrid > 0) {
                for (int ic = 1; ic <= constants::twotondim; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    int ixyz[3]; int idx = ic - 1; ixyz[2] = idx / 4; idx %= 4; ixyz[1] = idx / 2; ixyz[0] = idx % 2;
                    real_t x = grid_.xg[0 * grid_.ngridmax + igrid - 1] + (ixyz[0] - 0.5) * dx;
                    real_t y = (NDIM < 2) ? 0.5 : (grid_.xg[1 * grid_.ngridmax + igrid - 1] + (ixyz[1] - 0.5) * dx);
                    real_t z = (NDIM < 3) ? 0.5 : (grid_.xg[2 * grid_.ngridmax + igrid - 1] + (ixyz[2] - 0.5) * dx);
                    apply_ana(idc, x, y, z);
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    } else {
        for (int il = 1; il <= grid_.nlevelmax; ++il) {
            region_condinit(il);
        }
    }
}

void Initializer::load_grafic() {
    std::string path = config_.get("init_params", "initfile", "");
    if (path.empty()) path = config_.get("init_params", "initfile(1)", "");
    if (path.empty()) return;

    int lmin = params::levelmin;
    
    // 1. Read header from ic_velcx
    std::string header_file = path + "/ic_velcx";
    std::ifstream hfile(header_file, std::ios::binary);
    if (!hfile.is_open()) {
        std::cerr << "[Initializer] Could not open " << header_file << std::endl;
        return;
    }

    int32_t s; hfile.read((char*)&s, 4);
    int32_t n1, n2, n3; hfile.read((char*)&n1, 4); hfile.read((char*)&n2, 4); hfile.read((char*)&n3, 4);
    float dxini, xoff1, xoff2, xoff3, astart, om_m, om_l, h0_val;
    hfile.read((char*)&dxini, 4);
    hfile.read((char*)&xoff1, 4); hfile.read((char*)&xoff2, 4); hfile.read((char*)&xoff3, 4);
    hfile.read((char*)&astart, 4);
    hfile.read((char*)&om_m, 4); hfile.read((char*)&om_l, 4);
    hfile.read((char*)&h0_val, 4);
    hfile.read((char*)&s, 4);
    hfile.close();

    int npart_l = n1 * n2 * n3;
    std::cout << "[Initializer] Grafic Header: n=(" << n1 << "," << n2 << "," << n3 << ") astart=" << astart << std::endl;
    std::cout << "[Initializer] Loading " << npart_l << " particles (Level " << lmin << ")." << std::endl;

    grid_.resize_particles(npart_l);
    grid_.npart = npart_l;

    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (lmin - 1)));
    
    for (int i3 = 0; i3 < n3; ++i3) {
        for (int i2 = 0; i2 < n2; ++i2) {
            for (int i1 = 0; i1 < n1; ++i1) {
                int ip = i3 * n1 * n2 + i2 * n1 + i1;
                grid_.xp[0 * grid_.npartmax + ip] = (i1 + 0.5) * dx;
                if (NDIM > 1) grid_.xp[1 * grid_.npartmax + ip] = (i2 + 0.5) * dx;
                if (NDIM > 2) grid_.xp[2 * grid_.npartmax + ip] = (i3 + 0.5) * dx;

                grid_.mp[ip] = std::pow(0.5, 3 * lmin) * (1.0 - 0.045/0.3);
                grid_.levelp[ip] = lmin;
                grid_.idp[ip] = ip + 1;
            }
        }
    }

    // 2. Read velocities
    std::vector<std::string> v_files = {"/ic_velcx", "/ic_velcy", "/ic_velcz"};
    for (int d = 0; d < NDIM; ++d) {
        std::ifstream file(path + v_files[d], std::ios::binary);
        if (!file.is_open()) continue;
        
        file.read((char*)&s, 4); file.seekg(s, std::ios::cur); file.read((char*)&s, 4); // Skip header

        std::vector<float> plane(n1 * n2);
        for (int i3 = 0; i3 < n3; ++i3) {
            file.read((char*)&s, 4); file.read((char*)plane.data(), n1 * n2 * 4); file.read((char*)&s, 4);
            for (int i2 = 0; i2 < n2; ++i2) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int ip = i3 * n1 * n2 + i2 * n1 + i1;
                    if (ip < grid_.npartmax) {
                        grid_.vp[d * grid_.npartmax + ip] = plane[i2 * n1 + i1];
                    }
                }
            }
        }
    }

    std::cout << "[Initializer] Starting displacement check." << std::endl;
    // 3. Optional: Read positions (displacements)
    std::vector<std::string> x_files = {"/ic_poscx", "/ic_poscy", "/ic_poscz"};
    bool has_pos = true;
    for (int d = 0; d < NDIM; ++d) {
        std::ifstream file(path + x_files[d], std::ios::binary);
        if (!file.is_open()) { 
            std::cout << "[Initializer] Displacement file " << x_files[d] << " not found, will use Zeldovich." << std::endl;
            has_pos = false; break; 
        }
        
        file.read((char*)&s, 4); file.seekg(s, std::ios::cur); file.read((char*)&s, 4); // Skip header
        std::vector<float> plane(n1 * n2);
        for (int i3 = 0; i3 < n3; ++i3) {
            file.read((char*)&s, 4); file.read((char*)plane.data(), n1 * n2 * 4); file.read((char*)&s, 4);
            for (int i2 = 0; i2 < n2; ++i2) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int ip = i3 * n1 * n2 + i2 * n1 + i1;
                    if (ip < grid_.npartmax) grid_.xp[d * grid_.npartmax + ip] += plane[i2 * n1 + i1] / params::boxlen;
                }
            }
        }
    }

    if (!has_pos) {
        // Zeldovich approximation: xp = xp + vp
        for (int d = 0; d < NDIM; ++d) {
            for (int ip = 0; ip < grid_.npart; ++ip) {
                grid_.xp[d * grid_.npartmax + ip] += grid_.vp[d * grid_.npartmax + ip];
            }
        }
    }
    
    // Final periodic wrap
    for (int ip = 0; ip < grid_.npart; ++ip) {
        for (int d = 0; d < NDIM; ++d) {
            if (grid_.xp[d * grid_.npartmax + ip] < 0) grid_.xp[d * grid_.npartmax + ip] += 1.0;
            if (grid_.xp[d * grid_.npartmax + ip] >= 1.0) grid_.xp[d * grid_.npartmax + ip] -= 1.0;
        }
    }
}

void Initializer::region_condinit(int ilevel) {
    int nreg = config_.get_int("init_params", "nregion", 0);
    if (nreg == 0) return;

    real_t gam = grid_.gamma;
    
    std::vector<std::string> reg_type(nreg);
    std::vector<double> x_c(nreg, 0.5), y_c(nreg, 0.5), z_c(nreg, 0.5);
    std::vector<double> lx(nreg, 1.0), ly(nreg, 1.0), lz(nreg, 1.0);
    std::vector<double> dr(nreg, 1.0), ur(nreg, 0.0), vr(nreg, 0.0), wr(nreg, 0.0), pr(nreg, 1.0);
    std::vector<double> exp_reg(nreg, 10.0);

    auto get_list = [&](const std::string& block, const std::string& key) -> std::vector<std::string> {
        std::string val = config_.get(block, key, "");
        std::vector<std::string> res;
        std::stringstream ss(val);
        std::string item;
        while (std::getline(ss, item, ',')) {
            item = config_.trim(item);
            if (!item.empty() && (item.front() == '\'' || item.front() == '"')) item = item.substr(1);
            if (!item.empty() && (item.back() == '\'' || item.back() == '"')) item.pop_back();
            res.push_back(item);
        }
        return res;
    };

    auto dr_list = get_list("init_params", "d_region");
    auto ur_list = get_list("init_params", "u_region");
    auto vr_list = get_list("init_params", "v_region");
    auto wr_list = get_list("init_params", "w_region");
    auto pr_list = get_list("init_params", "p_region");
    auto var_list = get_list("init_params", "var_region");
    auto exp_list = get_list("init_params", "exp_region");
    auto xc_list = get_list("init_params", "x_center");
    auto yc_list = get_list("init_params", "y_center");
    auto zc_list = get_list("init_params", "z_center");
    auto lx_list = get_list("init_params", "length_x");
    auto ly_list = get_list("init_params", "length_y");
    auto lz_list = get_list("init_params", "length_z");

    auto reg_type_list = get_list("init_params", "region_type");

    for(int i=1; i<=nreg; ++i) {
        reg_type[i-1] = (i <= (int)reg_type_list.size()) ? reg_type_list[i-1] : "square";
        if (i <= (int)dr_list.size()) dr[i-1] = std::stod(dr_list[i-1]);
        if (i <= (int)ur_list.size()) ur[i-1] = std::stod(ur_list[i-1]);
        if (i <= (int)vr_list.size()) vr[i-1] = std::stod(vr_list[i-1]);
        if (i <= (int)wr_list.size()) wr[i-1] = std::stod(wr_list[i-1]);
        if (i <= (int)pr_list.size()) pr[i-1] = std::stod(pr_list[i-1]);
        if (i <= (int)exp_list.size()) exp_reg[i-1] = std::stod(exp_list[i-1]);
        if (i <= (int)xc_list.size()) x_c[i-1] = std::stod(xc_list[i-1]);
        if (i <= (int)yc_list.size()) y_c[i-1] = std::stod(yc_list[i-1]);
        if (i <= (int)zc_list.size()) z_c[i-1] = std::stod(zc_list[i-1]);
        if (i <= (int)lx_list.size()) lx[i-1] = std::stod(lx_list[i-1]);
        if (i <= (int)ly_list.size()) ly[i-1] = std::stod(ly_list[i-1]);
        if (i <= (int)lz_list.size()) lz[i-1] = std::stod(lz_list[i-1]);
    }

    auto apply_to_cell = [&](int idc, real_t x, real_t y, real_t z) {
        real_t dr_j = 1.0 + 1e-6 * (real_t(rand())/RAND_MAX - 0.5);
        real_t ur_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);
        real_t vr_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);
        real_t wr_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);

        // Default to Region 1 values if available, otherwise 0
        if (nreg > 0) {
            grid_.uold(idc, 1) = dr[0] * dr_j;
            grid_.uold(idc, 2) = dr[0] * (ur[0] + ur_j);
            grid_.uold(idc, 3) = dr[0] * (vr[0] + vr_j);
            grid_.uold(idc, 4) = dr[0] * (wr[0] + wr_j);
            real_t e_kin = 0.5 * dr[0] * (std::pow(ur[0]+ur_j, 2) + std::pow(vr[0]+vr_j, 2) + std::pow(wr[0]+wr_j, 2));
            real_t e_int = pr[0] / (gam - 1.0);
            grid_.uold(idc, 5) = e_kin + e_int;
            
            int nener = config_.get_int("hydro_params", "nener", 0);
            int npassive = 0;
            if (!var_list.empty()) npassive = var_list.size() / nreg;
            
            for(int ip=1; ip<=npassive; ++ip) grid_.uold(idc, 5 + nener + ip) = dr[0] * std::stod(var_list[ip-1]) * dr_j;
        }

        for (int ir = 0; ir < nreg; ++ir) {
            bool match = false;
            if (reg_type[ir] == "square") {
                real_t xn = 2.0 * std::abs(x - x_c[ir]) / lx[ir];
                real_t yn = (NDIM > 1) ? 2.0 * std::abs(y - y_c[ir]) / ly[ir] : 0.0;
                real_t zn = (NDIM > 2) ? 2.0 * std::abs(z - z_c[ir]) / lz[ir] : 0.0;
                real_t r = 0.0;
                if (exp_reg[ir] < 10.0) {
                    r = std::pow(std::pow(xn, exp_reg[ir]) + std::pow(yn, exp_reg[ir]) + std::pow(zn, exp_reg[ir]), 1.0 / exp_reg[ir]);
                } else {
                    r = std::max({xn, yn, zn});
                }
                if (r < 1.0) match = true;
            }
            if (match) {
                grid_.uold(idc, 1) = dr[ir] * dr_j;
                grid_.uold(idc, 2) = dr[ir] * (ur[ir] + ur_j);
                grid_.uold(idc, 3) = dr[ir] * (vr[ir] + vr_j);
                grid_.uold(idc, 4) = dr[ir] * (wr[ir] + wr_j);
                real_t e_kin = 0.5 * dr[ir] * (std::pow(ur[ir]+ur_j, 2) + std::pow(vr[ir]+vr_j, 2) + std::pow(wr[ir]+wr_j, 2));
                real_t e_int = pr[ir] / (gam - 1.0);
                grid_.uold(idc, 5) = e_kin + e_int;
                
                int nener = config_.get_int("hydro_params", "nener", 0);
                int npassive = 0;
                if (!var_list.empty()) npassive = var_list.size() / nreg;

                for(int ie=1; ie<=nener; ++ie) {
                    std::stringstream ss; ss << "prad_region(" << ie << "," << ir+1 << ")";
                    real_t p_rad = config_.get_double("init_params", ss.str(), 0.0);
                    if (p_rad > 0) {
                        real_t e_rad = p_rad / (gam - 1.0);
                        grid_.uold(idc, 5 + ie) = e_rad * dr_j;
                        grid_.uold(idc, 5) += e_rad * dr_j;
                    }
                }
                for(int ip=1; ip<=npassive; ++ip) {
                    grid_.uold(idc, 5 + nener + ip) = dr[ir] * std::stod(var_list[ir * npassive + (ip-1)]) * dr_j;
                }
            }
        }
    };

    if (ilevel == 1) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            int idx = i - 1;
            int iz = idx / (grid_.nx * grid_.ny); idx %= (grid_.nx * grid_.ny);
            int iy = idx / grid_.nx;
            int ix = idx % grid_.nx;
            real_t dx_coarse = grid_.boxlen / std::max({grid_.nx, grid_.ny, grid_.nz});
            real_t x = (ix + 0.5) * dx_coarse, y = (iy + 0.5) * dx_coarse, z = (iz + 0.5) * dx_coarse;
            apply_to_cell(i, x, y, z);
        }
    } else {
        int myid = MpiManager::instance().rank() + 1;
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t dx = grid_.boxlen / (real_t)(grid_.nx * (1 << (ilevel - 1)));
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                int ixyz[3] = {ix, iy, iz};
                real_t x = grid_.xg[0 * grid_.ngridmax + (igrid - 1)] + (ixyz[0] - 0.5) * dx;
                real_t y = (NDIM < 2) ? 0.5 : (grid_.xg[1 * grid_.ngridmax + (igrid - 1)] + (ixyz[1] - 0.5) * dx);
                real_t z = (NDIM < 3) ? 0.5 : (grid_.xg[2 * grid_.ngridmax + (igrid - 1)] + (ixyz[2] - 0.5) * dx);
                apply_to_cell(idc, x, y, z);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

} // namespace ramses
