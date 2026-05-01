#include "ramses/Initializer.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

namespace ramses {

void Initializer::apply_all() {
    std::string condinit_kind = config_.get("init_params", "condinit_kind", "");
    std::cout << "[Initializer] condinit_kind=\"" << condinit_kind << "\"" << std::endl;
    if (condinit_kind == "ana_disk_potential" || condinit_kind == "'ana_disk_potential'") {
        ana_disk_potential_condinit();
    } else if (condinit_kind == "orzag_tang" || condinit_kind == "'orzag_tang'") {
        orzag_tang_condinit();
    } else {
        region_condinit();
    }
}

void Initializer::orzag_tang_condinit() {}
void Initializer::ana_disk_potential_condinit() { region_condinit(); }

void Initializer::region_condinit() {
    int nreg = config_.get_int("init_params", "nregion", 1);
    real_t gam = config_.get_double("hydro_params", "gamma", 1.4);
    int nener = config_.get_int("hydro_params", "nener", 0);
    
    auto parse = [&](const std::string& k) {
        std::vector<real_t> v; std::string s = config_.get("init_params", k, "");
        if (s.empty()) return v;
        std::replace(s.begin(), s.end(), 'd', 'e'); std::replace(s.begin(), s.end(), 'D', 'e'); std::replace(s.begin(), s.end(), ',', ' ');
        std::stringstream ss(s); real_t val; while (ss >> val) v.push_back(val); return v;
    };

    auto drs = parse("d_region"), prs = parse("p_region"), urs = parse("u_region"), vrs = parse("v_region"), wrs = parse("w_region");
    auto xcs = parse("x_center"), lxs = parse("length_x"), ycs = parse("y_center"), lys = parse("length_y"), zcs = parse("z_center"), lzs = parse("length_z");

    real_t db = !drs.empty() ? drs[0] : 1.0, pb = !prs.empty() ? prs[0] : 1e-5, ub = !urs.empty() ? urs[0] : 0.0, vb = !vrs.empty() ? vrs[0] : 0.0, wb = !wrs.empty() ? wrs[0] : 0.0;
    int total_matches[10] = {0};

    auto apply_cell = [&](int id, int il) {
        real_t x_cell[3]; grid_.get_cell_center(id, x_cell);
        real_t xp = x_cell[0] * params::boxlen, yp = x_cell[1] * params::boxlen, zp = x_cell[2] * params::boxlen;
        
        real_t q[5] = {db, ub, vb, wb, pb};
        for (int r = 0; r < nreg; ++r) {
            std::string rkey = "region_type(" + std::to_string(r + 1) + ")";
            std::string rtype = config_.get("init_params", rkey, "square");
            real_t xc = (r < (int)xcs.size()) ? xcs[r] : 0.0, yc = (r < (int)ycs.size()) ? ycs[r] : 0.0, zc = (r < (int)zcs.size()) ? zcs[r] : 0.0;

            if (rtype == "square" || rtype == "'square'") {
                real_t lx = (r < (int)lxs.size()) ? lxs[r] : 1e10, ly = (r < (int)lys.size()) ? lys[r] : 1e10, lz = (r < (int)lzs.size()) ? lzs[r] : 1e10;
                real_t en = config_.get_double("init_params", "exp_region(" + std::to_string(r + 1) + ")", 10.0);
                real_t xn = 2.0 * std::abs(xp - xc) / lx, yn = 2.0 * std::abs(yp - yc) / ly, zn = 2.0 * std::abs(zp - zc) / lz;
                real_t rad = (en < 10.0) ? std::pow(std::pow(xn, en) + std::pow(yn, en) + std::pow(zn, en), 1.0/en) : std::max({xn, yn, zn});
                if (rad < 1.0) {
                    q[0] = (r < (int)drs.size()) ? drs[r] : db; q[1] = (r < (int)urs.size()) ? urs[r] : ub;
                    q[2] = (r < (int)vrs.size()) ? vrs[r] : vb; q[3] = (r < (int)wrs.size()) ? wrs[r] : wb;
                    q[4] = (r < (int)prs.size()) ? prs[r] : pb;
                    if (grid_.son[id] == 0 && r < 10) total_matches[r]++;
                }
            } else if (rtype == "point" || rtype == "'point'") {
                real_t dx = params::boxlen / static_cast<real_t>(params::nx * (1 << (il - 1)));
                real_t vol = std::pow(dx, NDIM);
                real_t xn = std::max(0.0, 1.0 - std::abs(xp - xc) / dx);
                real_t yn = std::max(0.0, 1.0 - std::abs(yp - yc) / dx);
                real_t zn = std::max(0.0, 1.0 - std::abs(zp - zc) / dx);
                real_t w = xn; if (NDIM > 1) w *= yn; if (NDIM > 2) w *= zn;
                if (w > 0.0) {
                    q[0] += ((r < (int)drs.size()) ? drs[r] : 0.0) * w / vol;
                    q[1] += ((r < (int)urs.size()) ? urs[r] : 0.0) * w;
                    q[2] += ((r < (int)vrs.size()) ? vrs[r] : 0.0) * w;
                    q[3] += ((r < (int)wrs.size()) ? wrs[r] : 0.0) * w;
                    q[4] += ((r < (int)prs.size()) ? prs[r] : 0.0) * w / vol;
                    if (grid_.son[id] == 0 && r < 10) total_matches[r]++;
                }
            }
        }
        q[0] = std::max(q[0], 1e-5);
        grid_.uold(id, 1) = q[0]; grid_.uold(id, 2) = q[0] * q[1]; grid_.uold(id, 3) = q[0] * q[2]; grid_.uold(id, 4) = q[0] * q[3];
        grid_.uold(id, 5) = q[4] / (gam - 1.0) + 0.5 * q[0] * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
        for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(id, iv) = grid_.uold(id, iv);
        grid_.cpu_map[id] = 1;
    };

    for (int i = 1; i <= grid_.ncell; ++i) {
        grid_.uold(i, 1) = db; grid_.uold(i, 2) = db * ub; grid_.uold(i, 3) = db * vb; grid_.uold(i, 4) = db * wb;
        grid_.uold(i, 5) = pb / (gam - 1.0) + 0.5 * db * (ub*ub + vb*vb + wb*wb);
        for (int iv = 1; iv <= grid_.nvar; ++iv) grid_.unew(i, iv) = grid_.uold(i, iv);
    }

    for (int i = 1; i <= grid_.ncoarse; ++i) apply_cell(i, 1);

    int myid = MpiManager::instance().rank() + 1;
    for (int il = 1; il <= grid_.nlevelmax; ++il) {
        int ig = grid_.get_headl(myid, il);
        while (ig > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                apply_cell(id, il);
            }
            ig = grid_.next[ig - 1];
        }
    }

    for(int r=0; r<nreg; ++r) if(total_matches[r] > 0) std::cout << "[Initializer] Region " << r+1 << " matched " << total_matches[r] << " cells." << std::endl;
}

void Initializer::initialize_cell(int i, int il, int nreg, real_t gam, int nener,
                                 const std::vector<real_t>& drs, const std::vector<real_t>& prs,
                                 const std::vector<real_t>& urs, const std::vector<real_t>& vrs,
                                 const std::vector<real_t>& wrs, const std::vector<real_t>& Ars,
                                 const std::vector<real_t>& Brs, const std::vector<real_t>& Crs,
                                 const std::vector<real_t>& xcs, const std::vector<real_t>& lxs,
                                 const std::vector<real_t>& ycs, const std::vector<real_t>& lys,
                                 const std::vector<real_t>& zcs, const std::vector<real_t>& lzs,
                                 const std::vector<std::vector<real_t>>& prads,
                                 real_t db, real_t pb, real_t ub, real_t vb, real_t wb,
                                 real_t Ab, real_t Bb, real_t Cb, int total_matches[]) {}

} // namespace ramses
