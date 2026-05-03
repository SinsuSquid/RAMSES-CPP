#include "ramses/Initializer.hpp"
#include "ramses/AmrGrid.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

namespace ramses {

void Initializer::apply_all() {
    for (int il = 1; il <= grid_.nlevelmax; ++il) {
        region_condinit(il);
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

    auto get_list = [&](const std::string& block, const std::string& key) -> std::vector<std::string> {
        std::string val = config_.get(block, key, "");
        std::vector<std::string> res;
        std::stringstream ss(val);
        std::string item;
        while (std::getline(ss, item, ',')) {
            res.push_back(config_.trim(item));
        }
        return res;
    };

    auto dr_list = get_list("init_params", "d_region");
    auto ur_list = get_list("init_params", "u_region");
    auto vr_list = get_list("init_params", "v_region");
    auto wr_list = get_list("init_params", "w_region");
    auto pr_list = get_list("init_params", "p_region");
    auto xc_list = get_list("init_params", "x_center");
    auto lx_list = get_list("init_params", "length_x");

    auto reg_type_list = get_list("init_params", "region_type");

    for(int i=1; i<=nreg; ++i) {
        reg_type[i-1] = (i <= (int)reg_type_list.size()) ? reg_type_list[i-1] : "square";
        if (i <= (int)dr_list.size()) dr[i-1] = std::stod(dr_list[i-1]);
        if (i <= (int)ur_list.size()) ur[i-1] = std::stod(ur_list[i-1]);
        if (i <= (int)vr_list.size()) vr[i-1] = std::stod(vr_list[i-1]);
        if (i <= (int)wr_list.size()) wr[i-1] = std::stod(wr_list[i-1]);
        if (i <= (int)pr_list.size()) pr[i-1] = std::stod(pr_list[i-1]);
        if (i <= (int)xc_list.size()) x_c[i-1] = std::stod(xc_list[i-1]);
        if (i <= (int)lx_list.size()) lx[i-1] = std::stod(lx_list[i-1]);
    }




    auto apply_to_cell = [&](int idc, real_t x, real_t y, real_t z) {
        for (int ir = 0; ir < nreg; ++ir) {
            bool match = false;
            if (reg_type[ir] == "square") {
                bool xm = std::abs(x - x_c[ir]) <= 0.5 * lx[ir] + 1e-10;
                bool ym = (NDIM < 2) || (std::abs(y - y_c[ir]) <= 0.5 * ly[ir] + 1e-10);
                bool zm = (NDIM < 3) || (std::abs(z - z_c[ir]) <= 0.5 * lz[ir] + 1e-10);
                if (xm && ym && zm) match = true;
            }
            if (match) {
                grid_.uold(idc, 1) = dr[ir];
                grid_.uold(idc, 2) = dr[ir] * ur[ir];
                grid_.uold(idc, 3) = dr[ir] * vr[ir];
                grid_.uold(idc, 4) = dr[ir] * wr[ir];
                real_t e_kin = 0.5 * dr[ir] * (ur[ir]*ur[ir] + vr[ir]*vr[ir] + wr[ir]*wr[ir]);
                real_t e_int = pr[ir] / (gam - 1.0);
                grid_.uold(idc, 5) = e_kin + e_int;
#ifdef RAMSES_NENER
                int nener = RAMSES_NENER;
#else
                int nener = config_.get_int("hydro_params", "nener", 0);
#endif
                for(int ie=1; ie<=nener; ++ie) {
                    std::stringstream ss; ss << "prad_region(" << ie << "," << ir+1 << ")";
                    real_t p_rad = config_.get_double("init_params", ss.str(), 0.0);
                    if (p_rad > 0) {
                        real_t e_rad = p_rad / (gam - 1.0);
                        grid_.uold(idc, 5 + ie) = e_rad;
                        grid_.uold(idc, 5) += e_rad;
                    }
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
        int myid = 1;
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
