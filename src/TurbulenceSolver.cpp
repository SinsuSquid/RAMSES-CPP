#include "ramses/TurbulenceSolver.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <random>

namespace ramses {

TurbulenceSolver::TurbulenceSolver(AmrGrid& grid, Config& config) 
    : grid_(grid), config_(config) {
    turb_gs_ = config_.get_int("turb_params", "turb_gs", 64);
    turb_min_rho_ = config_.get_double("turb_params", "turb_min_rho", 0.0);
    sol_frac_ = config_.get_double("turb_params", "sol_frac", 1.0);
}

TurbulenceSolver::~TurbulenceSolver() {}

void TurbulenceSolver::init() {
    size_t size = (size_t)NDIM * turb_gs_ * turb_gs_ * (NDIM == 3 ? turb_gs_ : 1);
    afield_now_.assign(size, 0.0);
    generate_forcing_field();
}

void TurbulenceSolver::update_fields(real_t dt) {
    // In a full implementation, this would update the OU process.
    // For now, we regenerate or maintain the field.
    // generate_forcing_field(); 
}

void TurbulenceSolver::apply_forcing(int ilevel, real_t dt) {
    int myid = MpiManager::instance().rank() + 1;
    int n2d_val = (1 << NDIM);
    int iener = NDIM + 2; // Index for total energy

    auto process_cell = [&](int idc) {
        if (grid_.son.at(idc - 1) > 0) return;
        real_t rho = grid_.uold(idc, 1);
        if (rho < turb_min_rho_) return;

        real_t xc[3];
        grid_.get_cell_center(idc, xc);

        real_t force[3] = {0};
        interpolate_force(xc, force);

        for (int idim = 0; idim < NDIM; ++idim) {
            real_t f = force[idim];
            // Update momentum: rho*v += rho*f*dt
            grid_.unew(idc, 2 + idim) += rho * f * dt;
            // Update energy: E += rho*v*f*dt + 0.5*rho*f^2*dt^2
            grid_.unew(idc, iener) += grid_.uold(idc, 2 + idim) * f * dt + 0.5 * rho * f * f * dt * dt;
        }
    };

    if (ilevel == 0) {
        for (int idc = 1; idc <= grid_.ncoarse; ++idc) process_cell(idc);
    }

    int igrid = grid_.get_headl(myid, ilevel);
    while (igrid > 0) {
        for (int ic = 1; ic <= n2d_val; ++ic) {
            int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
            process_cell(idc);
        }
        igrid = grid_.next[igrid - 1];
    }
}

void TurbulenceSolver::interpolate_force(const real_t x[3], real_t force[3]) {
    // Map x [0, boxlen] to [0, turb_gs]
    real_t r[3];
    for (int i = 0; i < NDIM; ++i) {
        r[i] = (x[i] / grid_.boxlen) * turb_gs_;
        while (r[i] < 0) r[i] += turb_gs_;
        while (r[i] >= turb_gs_) r[i] -= turb_gs_;
    }

    int bmin[3], bmax[3];
    real_t dr1[3], dr2[3];
    for (int i = 0; i < NDIM; ++i) {
        bmin[i] = (int)std::floor(r[i]);
        bmax[i] = (bmin[i] + 1) % turb_gs_;
        dr1[i] = r[i] - bmin[i];
        dr2[i] = 1.0 - dr1[i];
    }

    auto get_field = [&](int ix, int iy, int iz, int idim) {
        size_t idx;
        if (NDIM == 1) idx = (size_t)ix;
        else if (NDIM == 2) idx = (size_t)ix * turb_gs_ + iy;
        else idx = ((size_t)ix * turb_gs_ + iy) * turb_gs_ + iz;
        return afield_now_[idx * NDIM + idim];
    };

#if NDIM == 1
    for (int d = 0; d < NDIM; ++d) {
        force[d] = dr2[0] * get_field(bmin[0], 0, 0, d) +
                   dr1[0] * get_field(bmax[0], 0, 0, d);
    }
#elif NDIM == 2
    for (int d = 0; d < NDIM; ++d) {
        force[d] = dr2[0] * dr2[1] * get_field(bmin[0], bmin[1], 0, d) +
                   dr1[0] * dr2[1] * get_field(bmax[0], bmin[1], 0, d) +
                   dr2[0] * dr1[1] * get_field(bmin[0], bmax[1], 0, d) +
                   dr1[0] * dr1[1] * get_field(bmax[0], bmax[1], 0, d);
    }
#else
    for (int d = 0; d < NDIM; ++d) {
        force[d] = dr2[0] * dr2[1] * dr2[2] * get_field(bmin[0], bmin[1], bmin[2], d) +
                   dr1[0] * dr2[1] * dr2[2] * get_field(bmax[0], bmin[1], bmin[2], d) +
                   dr2[0] * dr1[1] * dr2[2] * get_field(bmin[0], bmax[1], bmin[2], d) +
                   dr1[0] * dr1[1] * dr2[2] * get_field(bmax[0], bmax[1], bmin[2], d) +
                   dr2[0] * dr2[1] * dr1[2] * get_field(bmin[0], bmin[1], bmax[2], d) +
                   dr1[0] * dr2[1] * dr1[2] * get_field(bmax[0], bmin[1], bmax[2], d) +
                   dr2[0] * dr1[1] * dr1[2] * get_field(bmin[0], bmax[1], bmax[2], d) +
                   dr1[0] * dr1[1] * dr1[2] * get_field(bmax[0], bmax[1], bmax[2], d);
    }
#endif
}

void TurbulenceSolver::generate_forcing_field() {
    // Mode-Sum Fallback Implementation
    // Generates a random solenoidal field using a sum of 100 Fourier modes
    int nmodes = 100;
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::uniform_real_distribution<real_t> dist(-1.0, 1.0);
    std::uniform_real_distribution<real_t> angle_dist(0.0, 2.0 * M_PI);

    struct Mode {
        real_t k[3];
        real_t A[3];
        real_t phi;
    };
    std::vector<Mode> modes(nmodes);

    for (int i = 0; i < nmodes; ++i) {
        // Random wavevector in [1, 4] range
        real_t kmag = 0;
        while (kmag < 1.0 || kmag > 4.0) {
            for (int d = 0; d < 3; ++d) modes[i].k[d] = std::floor(dist(gen) * 5.0);
            kmag = std::sqrt(modes[i].k[0]*modes[i].k[0] + modes[i].k[1]*modes[i].k[1] + modes[i].k[2]*modes[i].k[2]);
        }
        
        // Random amplitude vector
        real_t amp[3];
        for (int d = 0; d < 3; ++d) amp[d] = dist(gen);
        
        // Project to solenoidal (dot(k, A) = 0)
        real_t k_dot_amp = 0;
        for (int d = 0; d < 3; ++d) k_dot_amp += modes[i].k[d] * amp[d];
        for (int d = 0; d < 3; ++d) modes[i].A[d] = amp[d] - (k_dot_amp / (kmag*kmag)) * modes[i].k[d];
        
        modes[i].phi = angle_dist(gen);
    }

    // Fill the grid
    for (int ix = 0; ix < turb_gs_; ++ix) {
        for (int iy = 0; iy < (NDIM > 1 ? turb_gs_ : 1); ++iy) {
            for (int iz = 0; iz < (NDIM > 2 ? turb_gs_ : 1); ++iz) {
                real_t x[3] = {(real_t)ix / turb_gs_, (real_t)iy / turb_gs_, (real_t)iz / turb_gs_};
                size_t idx;
                if (NDIM == 1) idx = (size_t)ix;
                else if (NDIM == 2) idx = (size_t)ix * turb_gs_ + iy;
                else idx = ((size_t)ix * turb_gs_ + iy) * turb_gs_ + iz;

                for (int d = 0; d < NDIM; ++d) {
                    real_t val = 0;
                    for (int m = 0; m < nmodes; ++m) {
                        real_t k_dot_x = modes[m].k[0] * x[0] + modes[m].k[1] * x[1] + modes[m].k[2] * x[2];
                        val += modes[m].A[d] * std::cos(2.0 * M_PI * k_dot_x + modes[m].phi);
                    }
                    afield_now_[idx * NDIM + d] = val * 0.01; // Scale factor
                }
            }
        }
    }
}

} // namespace ramses
