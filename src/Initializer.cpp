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
    } else {
        for (int il = 1; il <= params::nlevelmax; ++il) {
            region_condinit(il);
        }
    }
}

void Initializer::region_condinit(int ilevel) {
    int nreg = config_.get_int("init_params", "nregion", 0);
    if (nreg == 0) return;

    real_t gam = grid_.gamma;
    std::vector<std::string> reg_type(nreg), d_list, p_list, ux_list, uy_list, uz_list, lx_list, ly_list, lz_list, exp_list, var_list;
    
    std::string rt_s = config_.get("init_params", "region_type", "");
    std::stringstream ssrt(rt_s); std::string item;
    for(int i=0; i<nreg && (ssrt >> item); ++i) reg_type[i] = item;
    
    auto split = [&](const std::string& key) {
        std::string s = config_.get("init_params", key, "");
        std::replace(s.begin(), s.end(), ',', ' ');
        std::stringstream ss(s); std::vector<std::string> res;
        std::string tmp; while(ss >> tmp) res.push_back(tmp);
        return res;
    };

    d_list = split("d_region"); p_list = split("p_region");
    ux_list = split("u_region"); uy_list = split("v_region"); uz_list = split("w_region");
    lx_list = split("length_x"); ly_list = split("length_y"); lz_list = split("length_z");
    exp_list = split("exp_region");
    var_list = split("var_region");

    std::vector<real_t> x_c(nreg, 0.5), y_c(nreg, 0.5), z_c(nreg, 0.5);
    std::vector<real_t> dr(nreg, 1.0), pr(nreg, 1.0), ur(nreg, 0.0), vr(nreg, 0.0), wr(nreg, 0.0);
    std::vector<real_t> lx(nreg, 1.0), ly(nreg, 1.0), lz(nreg, 1.0), exp_reg(nreg, 2.0);

    std::vector<std::string> xc_list = split("x_center"), yc_list = split("y_center"), zc_list = split("z_center");

    for (int i = 1; i <= nreg; ++i) {
        if (i <= (int)xc_list.size()) x_c[i-1] = std::stod(xc_list[i-1]);
        if (i <= (int)yc_list.size()) y_c[i-1] = std::stod(yc_list[i-1]);
        if (i <= (int)zc_list.size()) z_c[i-1] = std::stod(zc_list[i-1]);
        if (i <= (int)d_list.size()) dr[i-1] = std::stod(d_list[i-1]);
        if (i <= (int)p_list.size()) pr[i-1] = std::stod(p_list[i-1]);
        if (i <= (int)ux_list.size()) ur[i-1] = std::stod(ux_list[i-1]);
        if (i <= (int)uy_list.size()) vr[i-1] = std::stod(uy_list[i-1]);
        if (i <= (int)uz_list.size()) wr[i-1] = std::stod(uz_list[i-1]);
        if (i <= (int)exp_list.size()) exp_reg[i-1] = std::stod(exp_list[i-1]);
        if (i <= (int)lx_list.size()) lx[i-1] = std::stod(lx_list[i-1]);
        if (i <= (int)ly_list.size()) ly[i-1] = std::stod(ly_list[i-1]);
        if (i <= (int)lz_list.size()) lz[i-1] = std::stod(lz_list[i-1]);
    }

    auto apply_to_cell = [&](int idc, real_t x, real_t y, real_t z) {
        real_t dr_j = 1.0 + 1e-6 * (real_t(rand())/RAND_MAX - 0.5);
        real_t ur_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);
        real_t vr_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);
        real_t wr_j = 1e-2 * (real_t(rand())/RAND_MAX - 0.5);

        if (nreg > 0) {
            grid_.uold(idc, 1) = dr[0] * dr_j;
            grid_.uold(idc, 2) = dr[0] * (ur[0] + ur_j);
            if (NDIM > 1) grid_.uold(idc, 3) = dr[0] * (vr[0] + vr_j);
            if (NDIM > 2) grid_.uold(idc, 4) = dr[0] * (wr[0] + wr_j);
            
            real_t v2 = std::pow(ur[0]+ur_j, 2);
            if (NDIM > 1) v2 += std::pow(vr[0]+vr_j, 2);
            if (NDIM > 2) v2 += std::pow(wr[0]+wr_j, 2);
            real_t e_kin = 0.5 * dr[0] * v2;
            real_t e_int = pr[0] / (gam - 1.0);
            grid_.uold(idc, NDIM + 2) = e_kin + e_int;
            
            int nener = config_.get_int("hydro_params", "nener", 0);
            int npassive = 0;
            if (!var_list.empty()) npassive = var_list.size() / nreg;
            
            for(int ip=1; ip<=npassive; ++ip) grid_.uold(idc, NDIM + 2 + nener + ip) = dr[0] * std::stod(var_list[ip-1]) * dr_j;
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
                if (NDIM > 1) grid_.uold(idc, 3) = dr[ir] * (vr[ir] + vr_j);
                if (NDIM > 2) grid_.uold(idc, 4) = dr[ir] * (wr[ir] + wr_j);
                
                real_t v2 = std::pow(ur[ir]+ur_j, 2);
                if (NDIM > 1) v2 += std::pow(vr[ir]+vr_j, 2);
                if (NDIM > 2) v2 += std::pow(wr[ir]+wr_j, 2);
                real_t e_kin = 0.5 * dr[ir] * v2;
                real_t e_int = pr[ir] / (gam - 1.0);
                grid_.uold(idc, NDIM + 2) = e_kin + e_int;
                
                int nener = config_.get_int("hydro_params", "nener", 0);
                int npassive = 0;
                if (!var_list.empty()) npassive = var_list.size() / nreg;

                for(int ie=1; ie<=nener; ++ie) {
                    std::stringstream ss; ss << "prad_region(" << ie << "," << ir+1 << ")";
                    real_t p_rad = config_.get_double("init_params", ss.str(), 0.0);
                    if (p_rad > 0) {
                        real_t e_rad = p_rad / (gam - 1.0);
                        grid_.uold(idc, NDIM + 2 + ie) = e_rad * dr_j;
                        grid_.uold(idc, NDIM + 2) += e_rad * dr_j;
                    }
                }
                for(int ip=1; ip<=npassive; ++ip) {
                    grid_.uold(idc, NDIM + 2 + nener + ip) = dr[ir] * std::stod(var_list[ir * npassive + (ip-1)]) * dr_j;
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

void Initializer::load_grafic() {
    std::cout << "[Initializer] Loading Grafic ICs..." << std::endl;
}

} // namespace ramses
