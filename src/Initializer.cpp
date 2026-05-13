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
#include <random>

namespace ramses {

Initializer::~Initializer() {}

void Initializer::apply_all() {
    namespace p = ramses::params;
    std::string filetype = config_.get("init_params", "filetype", "region");
    if (filetype == "grafic") {
        load_grafic();
    } else {
        for (int il = 0; il <= params::nlevelmax; ++il) {
            region_condinit(il);
        }
    }
}

void Initializer::region_condinit(int ilevel) {
    int nreg = config_.get_int("init_params", "nregion", 0);
    if (nreg == 0) return;

    real_t gam = grid_.gamma;
    
    auto split = [&](const std::string& key) {
        std::string s = config_.get("init_params", key, "");
        std::replace(s.begin(), s.end(), ',', ' ');
        std::stringstream ss(s); std::vector<std::string> res;
        std::string tmp; while(ss >> tmp) res.push_back(tmp);
        return res;
    };

    std::vector<std::string> rt_list = split("region_type");
    std::vector<std::string> d_list = split("d_region");
    std::vector<std::string> p_list = split("p_region");
    std::vector<std::string> ux_list = split("u_region");
    std::vector<std::string> uy_list = split("v_region");
    std::vector<std::string> uz_list = split("w_region");
    std::vector<std::string> lx_list = split("length_x");
    std::vector<std::string> ly_list = split("length_y");
    std::vector<std::string> lz_list = split("length_z");
    std::vector<std::string> xc_list = split("x_center");
    std::vector<std::string> yc_list = split("y_center");
    std::vector<std::string> zc_list = split("z_center");
    std::vector<std::string> exp_list = split("exp_region");
    std::vector<std::string> var_list = split("var_region");

    std::vector<real_t> x_c(nreg, 0.5), y_c(nreg, 0.5), z_c(nreg, 0.5);
    std::vector<real_t> dr(nreg, 1.0), pr(nreg, 1.0), ur(nreg, 0.0), vr(nreg, 0.0), wr(nreg, 0.0);
    std::vector<real_t> lx(nreg, 1.0), ly(nreg, 1.0), lz(nreg, 1.0), exp_reg(nreg, 2.0);
    std::vector<std::string> reg_type(nreg, "cube");

    for (int i = 0; i < nreg; ++i) {
        if (i < (int)rt_list.size()) reg_type[i] = rt_list[i];
        if (i < (int)xc_list.size()) x_c[i] = std::stod(xc_list[i]);
        if (i < (int)yc_list.size()) y_c[i] = std::stod(yc_list[i]);
        if (i < (int)zc_list.size()) z_c[i] = std::stod(zc_list[i]);
        if (i < (int)d_list.size()) dr[i] = std::stod(d_list[i]);
        if (i < (int)p_list.size()) pr[i] = std::stod(p_list[i]);
        if (i < (int)ux_list.size()) ur[i] = std::stod(ux_list[i]);
        if (i < (int)uy_list.size()) vr[i] = std::stod(uy_list[i]);
        if (i < (int)uz_list.size()) wr[i] = std::stod(uz_list[i]);
        if (i < (int)lx_list.size()) lx[i] = std::stod(lx_list[i]);
        if (i < (int)ly_list.size()) ly[i] = std::stod(ly_list[i]);
        if (i < (int)lz_list.size()) lz[i] = std::stod(lz_list[i]);
        if (i < (int)exp_list.size()) exp_reg[i] = std::stod(exp_list[i]);
    }

    auto apply_to_cell = [&](int idc, real_t x, real_t y, real_t z) {
        for (int ir = 0; ir < nreg; ++ir) {
            bool in_reg = false;
            if (reg_type[ir] == "cube" || reg_type[ir] == "square") {
                if (std::abs(x - x_c[ir]) <= 0.5 * lx[ir] + 1e-10 &&
                    std::abs(y - y_c[ir]) <= 0.5 * ly[ir] + 1e-10 &&
                    std::abs(z - z_c[ir]) <= 0.5 * lz[ir] + 1e-10) in_reg = true;
            } else if (reg_type[ir] == "sphere") {
                real_t r2 = std::pow(x - x_c[ir], 2) + std::pow(y - y_c[ir], 2) + std::pow(z - z_c[ir], 2);
                if (r2 <= std::pow(lx[ir], 2) + 1e-10) in_reg = true;
            }

            if (in_reg) {
                grid_.uold(idc, 1) = dr[ir];
                grid_.uold(idc, 2) = dr[ir] * ur[ir];
                if (NDIM > 1) grid_.uold(idc, 3) = dr[ir] * vr[ir];
                if (NDIM > 2) grid_.uold(idc, 4) = dr[ir] * wr[ir];
                
                real_t v2 = ur[ir]*ur[ir];
                if (NDIM > 1) v2 += vr[ir]*vr[ir];
                if (NDIM > 2) v2 += wr[ir]*wr[ir];
                real_t e_kin = 0.5 * dr[ir] * v2;
                real_t e_int = pr[ir] / (gam - 1.0);
                grid_.uold(idc, NDIM + 2) = e_kin + e_int;
                
                int nener = config_.get_int("hydro_params", "nener", 0);
                for(int ip=1; ip<=grid_.nvar - (NDIM + 2 + nener); ++ip) {
                    grid_.uold(idc, NDIM + 2 + nener + ip) = 0.0;
                }
            }
        }
    };

    if (ilevel == 0) {
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            real_t xc[3]; grid_.get_cell_center(i, xc);
            apply_to_cell(i, xc[0], xc[1], xc[2]);
        }
    } else {
        int myid = MpiManager::instance().rank() + 1;
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                real_t xc[3]; grid_.get_cell_center(idc, xc);
                apply_to_cell(idc, xc[0], xc[1], xc[2]);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void Initializer::load_grafic() {
    std::cout << "[Initializer] Loading Grafic ICs..." << std::endl;
}

void Initializer::init_tracers() {
    namespace p = ramses::params;
    if (!p::tracer) return;

    if (p::tracer_feed_fmt == "inplace") {
        if (MpiManager::instance().rank() == 0) {
            std::cout << "[Initializer] Creating tracers 'inplace' based on gas density..." << std::endl;
        }

        int myid = MpiManager::instance().rank() + 1;
        int n2d = (1 << NDIM);
        
        // Use a fixed seed for reproducible tracers
        std::mt19937 gen(1337 + myid);
        std::uniform_real_distribution<real_t> dis(0.0, 1.0);

        for (int ilevel = p::levelmin; ilevel <= p::nlevelmax; ++ilevel) {
            real_t dx = p::boxlen / (real_t)(p::nx * (1 << ilevel));
            real_t vol = std::pow(dx, NDIM);

            auto process_cell = [&](int idc) {
                if (grid_.son[idc - 1] != 0) return;

                real_t mass = grid_.uold(idc, 1) * vol;
                if (p::tracer_mass <= 0) return;
                real_t n_real = mass / p::tracer_mass;
                int n_int = static_cast<int>(n_real);
                if (dis(gen) < (n_real - n_int)) n_int++;

                if (n_int > 0) {
                    real_t xc[3];
                    grid_.get_cell_center(idc, xc);
                    for (int i = 0; i < n_int; ++i) {
                        int ip = grid_.get_free_particle();
                        if (ip == 0) {
                            return;
                        }
                        for (int idim = 0; idim < NDIM; ++idim) {
                            grid_.xp[idim * grid_.npartmax + ip - 1] = xc[idim];
                            grid_.vp[idim * grid_.npartmax + ip - 1] = 0.0;
                        }
                        grid_.mp[ip - 1] = p::tracer_mass;
                        grid_.levelp[ip - 1] = ilevel;
                        grid_.family[ip - 1] = FAM_TRACER;
                        grid_.tag[ip - 1] = 0;
                        grid_.idp[ip - 1] = ip; // Temporary ID
                    }
                }
            };

            if (ilevel == p::levelmin) {
                for (int i = 1; i <= grid_.ncoarse; ++i) process_cell(i);
            }

            int ig = grid_.get_headl(myid, ilevel);
            while (ig > 0) {
                for (int ic = 1; ic <= n2d; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                    process_cell(idc);
                }
                ig = grid_.next[ig - 1];
            }
        }
    } else {
        if (MpiManager::instance().rank() == 0) {
            std::cout << "[Initializer] Warning: tracer_feed_fmt '" << p::tracer_feed_fmt << "' not implemented yet." << std::endl;
        }
    }
}

} // namespace ramses
