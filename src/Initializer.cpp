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
    int myid = 1;
    int igrid = grid_.get_headl(myid, ilevel);

    std::vector<std::string> reg_type(nreg);
    std::vector<double> x_c(nreg, 0.5), y_c(nreg, 0.5), z_c(nreg, 0.5);
    std::vector<double> lx(nreg, 1.0), ly(nreg, 1.0), lz(nreg, 1.0);
    std::vector<double> dr(nreg, 1.0), ur(nreg, 0.0), vr(nreg, 0.0), wr(nreg, 0.0), pr(nreg, 1.0);

    for(int i=1; i<=nreg; ++i) {
        std::stringstream ss; ss << "(" << i << ")"; std::string idx = ss.str();
        reg_type[i-1] = config_.get("init_params", "region_type" + idx, "square");
        x_c[i-1] = config_.get_double("init_params", "x_center" + idx, 0.5);
        y_c[i-1] = config_.get_double("init_params", "y_center" + idx, 0.5);
        z_c[i-1] = config_.get_double("init_params", "z_center" + idx, 0.5);
        lx[i-1] = config_.get_double("init_params", "length_x" + idx, 1.0);
        ly[i-1] = config_.get_double("init_params", "length_y" + idx, 1.0);
        lz[i-1] = config_.get_double("init_params", "length_z" + idx, 1.0);
        dr[i-1] = config_.get_double("init_params", "d_region" + idx, 1.0);
        ur[i-1] = config_.get_double("init_params", "u_region" + idx, 0.0);
        vr[i-1] = config_.get_double("init_params", "v_region" + idx, 0.0);
        wr[i-1] = config_.get_double("init_params", "w_region" + idx, 0.0);
        pr[i-1] = config_.get_double("init_params", "p_region" + idx, 1.0);
    }

    while (igrid > 0) {
        for (int ic = 1; ic <= constants::twotondim; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            real_t dx = grid_.boxlen / static_cast<real_t>(1 << (ilevel - 1));
            int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
            real_t x = grid_.xg[0 * grid_.ngridmax + (igrid - 1)] + (ix - 0.5) * dx;
            real_t y = grid_.xg[1 * grid_.ngridmax + (igrid - 1)] + (iy - 0.5) * dx;
            real_t z = grid_.xg[2 * grid_.ngridmax + (igrid - 1)] + (iz - 0.5) * dx;

            for (int ir = 0; ir < nreg; ++ir) {
                bool match = false;
                if (reg_type[ir] == "square") {
                    if (std::abs(x - x_c[ir]) <= 0.5 * lx[ir] + 1e-10 &&
                        std::abs(y - y_c[ir]) <= 0.5 * ly[ir] + 1e-10 &&
                        std::abs(z - z_c[ir]) <= 0.5 * lz[ir] + 1e-10) match = true;
                }
                if (match) {
                    grid_.uold(idc, 1) = dr[ir];
                    grid_.uold(idc, 2) = dr[ir] * ur[ir];
                    grid_.uold(idc, 3) = dr[ir] * vr[ir];
                    grid_.uold(idc, 4) = dr[ir] * wr[ir];
                    real_t e_kin = 0.5 * dr[ir] * (ur[ir] * ur[ir] + vr[ir] * vr[ir] + wr[ir] * wr[ir]);
                    real_t e_int = pr[ir] / (gam - 1.0);
                    grid_.uold(idc, 5) = e_kin + e_int;
                    for(int ie=1; ie<=10; ++ie) {
                        std::stringstream ss; ss << "prad_region(" << ir+1 << "," << ie << ")";
                        real_t p_rad = config_.get_double("init_params", ss.str(), 0.0);
                        if (p_rad > 0) {
                            real_t e_rad = p_rad / (gam - 1.0);
                            grid_.uold(idc, 5 + ie) = e_rad;
                            grid_.uold(idc, 5) += e_rad;
                        }
                    }
                }
            }
        }
        igrid = grid_.next[igrid - 1];
    }
}

} // namespace ramses
