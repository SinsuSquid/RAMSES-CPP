#include "ramses/ParticleSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

ParticleSolver::~ParticleSolver() {}

void ParticleSolver::assign_mass(int ilevel) {
    assign_mass_fine(ilevel);
}

void ParticleSolver::assign_mass_fine(int ilevel) {
    if (grid_.npart == 0) return;

    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
    real_t vol = std::pow(dx, NDIM);
    real_t inv_dx = 1.0 / dx;
    
    int myid = MpiManager::instance().rank() + 1;
    for (int igrid = grid_.get_headl(myid, ilevel); igrid > 0; igrid = grid_.next[igrid - 1]) {
        for (int ip = grid_.headp[igrid - 1]; ip > 0; ip = grid_.nextp[ip - 1]) {
            if (grid_.family[ip - 1] == FAM_TRACER) continue;

            real_t xp = grid_.xp[0 * grid_.npartmax + ip - 1];
            real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;
            real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;

            int icell = find_cell_by_coords(xp, yp, zp);
            if (icell <= 0) continue;

            int ig = (icell > grid_.ncoarse) ? ((icell - grid_.ncoarse - 1) % grid_.ngridmax) + 1 : 0;
            int ic = (icell > grid_.ncoarse) ? ((icell - grid_.ncoarse - 1) / grid_.ngridmax) + 1 : icell;
            
            real_t xgc, ygc, zgc;
            if (ig > 0) {
                xgc = grid_.xg[0 * grid_.ngridmax + ig - 1];
                ygc = (NDIM > 1) ? grid_.xg[1 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
                zgc = (NDIM > 2) ? grid_.xg[2 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                xgc += (ix - 0.5) * dx; ygc += (iy - 0.5) * dx; zgc += (iz - 0.5) * dx;
            } else {
                real_t dx_coarse = params::boxlen / std::max({params::nx, params::ny, params::nz});
                int ixc = (ic - 1) % params::nx, iyc = ((ic - 1) / params::nx) % params::ny, izc = (ic - 1) / (params::nx * params::ny);
                xgc = (ixc + 0.5) * dx_coarse; ygc = (iyc + 0.5) * dx_coarse; zgc = (izc + 0.5) * dx_coarse;
            }

            real_t rx = (xp - xgc) * inv_dx;
            real_t ry = (NDIM > 1) ? (yp - ygc) * inv_dx : 0.0;
            real_t rz = (NDIM > 2) ? (zp - zgc) * inv_dx : 0.0;

            int nbors[27];
            grid_.get_27_cell_neighbors(icell, nbors);

            int sx = (rx > 0) ? 1 : -1;
            int sy = (ry > 0) ? 1 : -1;
            int sz = (rz > 0) ? 1 : -1;

            real_t wx1 = 1.0 - std::abs(rx), wx2 = std::abs(rx);
            real_t wy1 = 1.0 - std::abs(ry), wy2 = std::abs(ry);
            real_t wz1 = 1.0 - std::abs(rz), wz2 = std::abs(rz);

            auto assign = [&](int off_x, int off_y, int off_z, real_t w) {
                int idx = 13 + off_x + 3 * off_y + 9 * off_z;
                if (idx < 0 || idx >= 27) return;
                int target = nbors[idx];
                if (target > 0 && target <= (int)grid_.rho.size()) {
                    grid_.rho[target - 1] += grid_.mp[ip - 1] * w / vol;
                }
            };

            if (NDIM == 1) { assign(0, 0, 0, wx1); assign(sx, 0, 0, wx2); }
            else if (NDIM == 2) {
                assign(0, 0, 0, wx1 * wy1); assign(sx, 0, 0, wx2 * wy1);
                assign(0, sy, 0, wx1 * wy2); assign(sx, sy, 0, wx2 * wy2);
            } else {
                assign(0, 0, 0, wx1 * wy1 * wz1); assign(sx, 0, 0, wx2 * wy1 * wz1);
                assign(0, sy, 0, wx1 * wy2 * wz1); assign(sx, sy, 0, wx2 * wy2 * wz1);
                assign(0, 0, sz, wx1 * wy1 * wz2); assign(sx, 0, sz, wx2 * wy1 * wz2);
                assign(0, sy, sz, wx1 * wy2 * wz2); assign(sx, sy, sz, wx2 * wy2 * wz2);
            }
        }
    }
}

void ParticleSolver::move_fine(int ilevel, real_t dt) {
    if (grid_.npart == 0) return;

    real_t dx = params::boxlen / (real_t)(params::nx * (1 << (ilevel - 1)));
    real_t inv_dx = 1.0 / dx;
    int myid = MpiManager::instance().rank() + 1;

    for (int igrid = grid_.get_headl(myid, ilevel); igrid > 0; igrid = grid_.next[igrid - 1]) {
        for (int ip = grid_.headp[igrid - 1]; ip > 0; ip = grid_.nextp[ip - 1]) {
            real_t xp = grid_.xp[0 * grid_.npartmax + ip - 1];
            real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;
            real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;

            int icell = find_cell_by_coords(xp, yp, zp);
            if (icell <= 0) continue;

            int ig = (icell > grid_.ncoarse) ? ((icell - grid_.ncoarse - 1) % grid_.ngridmax) + 1 : 0;
            int ic = (icell > grid_.ncoarse) ? ((icell - grid_.ncoarse - 1) / grid_.ngridmax) + 1 : icell;
            
            real_t xgc, ygc, zgc;
            if (ig > 0) {
                xgc = grid_.xg[0 * grid_.ngridmax + ig - 1];
                ygc = (NDIM > 1) ? grid_.xg[1 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
                zgc = (NDIM > 2) ? grid_.xg[2 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
                int ix = (ic - 1) & 1, iy = ((ic - 1) & 2) >> 1, iz = ((ic - 1) & 4) >> 2;
                xgc += (ix - 0.5) * dx; ygc += (iy - 0.5) * dx; zgc += (iz - 0.5) * dx;
            } else {
                real_t dx_coarse = params::boxlen / std::max({params::nx, params::ny, params::nz});
                int ixc = (ic - 1) % params::nx, iyc = ((ic - 1) / params::nx) % params::ny, izc = (ic - 1) / (params::nx * params::ny);
                xgc = (ixc + 0.5) * dx_coarse; ygc = (iyc + 0.5) * dx_coarse; zgc = (izc + 0.5) * dx_coarse;
            }

            real_t rx = (xp - xgc) * inv_dx;
            real_t ry = (NDIM > 1) ? (yp - ygc) * inv_dx : 0.0;
            real_t rz = (NDIM > 2) ? (zp - zgc) * inv_dx : 0.0;

            int nbors[27];
            grid_.get_27_cell_neighbors(icell, nbors);

            int sx = (rx > 0) ? 1 : -1;
            int sy = (ry > 0) ? 1 : -1;
            int sz = (rz > 0) ? 1 : -1;

            real_t wx1 = 1.0 - std::abs(rx), wx2 = std::abs(rx);
            real_t wy1 = 1.0 - std::abs(ry), wy2 = std::abs(ry);
            real_t wz1 = 1.0 - std::abs(rz), wz2 = std::abs(rz);

            real_t f_interp[3] = {0, 0, 0};
            auto interpolate = [&](int off_x, int off_y, int off_z, real_t w) {
                int idx = 13 + off_x + 3 * off_y + 9 * off_z;
                if (idx < 0 || idx >= 27) return;
                int target = nbors[idx];
                if (target > 0 && target <= (int)grid_.rho.size()) {
                    for (int idim = 1; idim <= NDIM; ++idim) {
                        f_interp[idim - 1] += grid_.f(target, idim) * w;
                    }
                }
            };

            if (NDIM == 1) { interpolate(0, 0, 0, wx1); interpolate(sx, 0, 0, wx2); }
            else if (NDIM == 2) {
                interpolate(0, 0, 0, wx1 * wy1); interpolate(sx, 0, 0, wx2 * wy1);
                interpolate(0, sy, 0, wx1 * wy2); interpolate(sx, sy, 0, wx2 * wy2);
            } else {
                interpolate(0, 0, 0, wx1 * wy1 * wz1); interpolate(sx, 0, 0, wx2 * wy1 * wz1);
                interpolate(0, sy, 0, wx1 * wy2 * wz1); interpolate(sx, sy, 0, wx2 * wy2 * wz1);
                interpolate(0, 0, sz, wx1 * wy1 * wz2); interpolate(sx, 0, sz, wx2 * wy1 * wz2);
                interpolate(0, sy, sz, wx1 * wy2 * wz2); interpolate(sx, sy, sz, wx2 * wy2 * wz2);
            }

            if (grid_.family[ip - 1] == FAM_TRACER) {
                // Classical tracers: follow the gas velocity
                // We need to interpolate gas velocity from the grid.
                // Currently f_interp contains the force (from grid_.f).
                // We need a way to get velocity.
                // In RAMSES, velocity is uold(2:4)/uold(1).
                
                real_t v_gas[3] = {0, 0, 0};
                auto interp_v = [&](int off_x, int off_y, int off_z, real_t w) {
                    int idx = 13 + off_x + 3 * off_y + 9 * off_z;
                    if (idx < 0 || idx >= 27) return;
                    int target = nbors[idx];
                    if (target > 0) {
                        real_t rho_val = std::max(grid_.uold(target, 1), 1e-10);
                        for (int idim = 0; idim < NDIM; ++idim) {
                            v_gas[idim] += (grid_.uold(target, 2 + idim) / rho_val) * w;
                        }
                    }
                };

                if (NDIM == 1) { interp_v(0, 0, 0, wx1); interp_v(sx, 0, 0, wx2); }
                else if (NDIM == 2) {
                    interp_v(0, 0, 0, wx1 * wy1); interp_v(sx, 0, 0, wx2 * wy1);
                    interp_v(0, sy, 0, wx1 * wy2); interp_v(sx, sy, 0, wx2 * wy2);
                } else {
                    interp_v(0, 0, 0, wx1 * wy1 * wz1); interp_v(sx, 0, 0, wx2 * wy1 * wz1);
                    interp_v(0, sy, 0, wx1 * wy2 * wz1); interp_v(sx, sy, 0, wx2 * wy2 * wz1);
                    interp_v(0, 0, sz, wx1 * wy1 * wz2); interp_v(sx, 0, sz, wx2 * wy1 * wz2);
                    interp_v(0, sy, sz, wx1 * wy2 * wz2); interp_v(sx, sy, sz, wx2 * wy2 * wz2);
                }

                for (int idim = 0; idim < NDIM; ++idim) {
                    grid_.vp[idim * grid_.npartmax + ip - 1] = v_gas[idim];
                    grid_.xp[idim * grid_.npartmax + ip - 1] += grid_.vp[idim * grid_.npartmax + ip - 1] * dt;
                }
            } else {
                for (int idim = 1; idim <= NDIM; ++idim) {
                    grid_.vp[(idim - 1) * grid_.npartmax + ip - 1] += f_interp[idim - 1] * dt;
                    grid_.xp[(idim - 1) * grid_.npartmax + ip - 1] += grid_.vp[(idim - 1) * grid_.npartmax + ip - 1] * dt;
                }
            }

            for (int idim = 0; idim < NDIM; ++idim) {
                if (grid_.xp[idim * grid_.npartmax + ip - 1] < 0) grid_.xp[idim * grid_.npartmax + ip - 1] += params::boxlen;
                if (grid_.xp[idim * grid_.npartmax + ip - 1] >= params::boxlen) grid_.xp[idim * grid_.npartmax + ip - 1] -= params::boxlen;
            }
        }
    }
}

void ParticleSolver::move_particles(const std::vector<int>& ind_part, real_t dt) {
    // Standard move logic
}

void ParticleSolver::exchange_particles() {
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) return;

    int myid = mpi.rank() + 1;
    std::map<int, std::vector<ParticlePacket>> send_queues;
    std::vector<int> particles_to_free;

    for (int ilevel = 1; ilevel <= grid_.nlevelmax; ilevel++) {
        int igrid = grid_.get_headl(myid, ilevel);
        while (igrid > 0) {
            int ip = grid_.headp[igrid - 1];
            while (ip > 0) {
                int next_p = grid_.nextp[ip - 1];
                
                real_t xp = grid_.xp[0 * grid_.npartmax + ip - 1];
                real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;
                real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;
                
                int icell = find_cell_by_coords(xp, yp, zp);
                int target_cpu = grid_.cpu_map[icell - 1];
                
                if (target_cpu != myid) {
                    ParticlePacket p;
                    p.xp[0] = xp; p.xp[1] = yp; p.xp[2] = zp;
                    p.vp[0] = grid_.vp[0 * grid_.npartmax + ip - 1];
                    p.vp[1] = (NDIM > 1) ? grid_.vp[1 * grid_.npartmax + ip - 1] : 0.0;
                    p.vp[2] = (NDIM > 2) ? grid_.vp[2 * grid_.npartmax + ip - 1] : 0.0;
                    p.mp = grid_.mp[ip - 1];
                    p.idp = grid_.idp[ip - 1];
                    p.levelp = grid_.levelp[ip - 1];
                    p.family = grid_.family[ip - 1];
                    p.tag = grid_.tag[ip - 1];
                    
                    send_queues[target_cpu].push_back(p);
                    particles_to_free.push_back(ip);
                    
                    int prev_idx = grid_.prevp[ip - 1];
                    int next_idx = grid_.nextp[ip - 1];
                    if (prev_idx > 0) grid_.nextp[prev_idx - 1] = next_idx;
                    else grid_.headp[igrid - 1] = next_idx;
                    if (next_idx > 0) grid_.prevp[next_idx - 1] = prev_idx;
                    else grid_.tailp[igrid - 1] = prev_idx;
                    grid_.numbp[igrid - 1]--;
                }
                ip = next_p;
            }
            igrid = grid_.next[igrid - 1];
        }
    }

    for (int ip : particles_to_free) grid_.free_particle(ip);

#ifdef RAMSES_USE_MPI
    int ncpu = mpi.size();
    std::vector<int> send_counts(ncpu, 0), recv_counts(ncpu, 0);
    for (auto const& [cpu, queue] : send_queues) send_counts[cpu - 1] = queue.size();

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<MPI_Request> requests;
    std::vector<std::vector<ParticlePacket>> recv_queues(ncpu);

    for (int i = 0; i < ncpu; ++i) {
        if (i == mpi.rank()) continue;
        if (send_counts[i] > 0) {
            MPI_Request req;
            MPI_Isend(send_queues[i + 1].data(), send_counts[i] * sizeof(ParticlePacket), MPI_BYTE, i, 1, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
        if (recv_counts[i] > 0) {
            recv_queues[i].resize(recv_counts[i]);
            MPI_Request req;
            MPI_Irecv(recv_queues[i].data(), recv_counts[i] * sizeof(ParticlePacket), MPI_BYTE, i, 1, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
    }
    if (!requests.empty()) MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    for (int i = 0; i < ncpu; ++i) {
        for (const auto& p : recv_queues[i]) {
            int ip = grid_.get_free_particle();
            if (ip == 0) {
                std::cerr << "[ParticleSolver] Error: Out of particle slots during MPI exchange!" << std::endl;
                continue;
            }
            grid_.xp[0 * grid_.npartmax + ip - 1] = p.xp[0];
            grid_.xp[1 * grid_.npartmax + ip - 1] = p.xp[1];
            grid_.xp[2 * grid_.npartmax + ip - 1] = p.xp[2];
            grid_.vp[0 * grid_.npartmax + ip - 1] = p.vp[0];
            grid_.vp[1 * grid_.npartmax + ip - 1] = p.vp[1];
            grid_.vp[2 * grid_.npartmax + ip - 1] = p.vp[2];
            grid_.mp[ip - 1] = p.mp;
            grid_.idp[ip - 1] = p.idp;
            grid_.levelp[ip - 1] = p.levelp;
            grid_.family[ip - 1] = p.family;
            grid_.tag[ip - 1] = p.tag;
        }
    }
#endif

    relink();
}

void ParticleSolver::relink() {
    if (grid_.npart == 0) return;

    for (int ig = 1; ig <= grid_.ngridmax; ++ig) {
        grid_.headp[ig - 1] = 0;
        grid_.tailp[ig - 1] = 0;
        grid_.numbp[ig - 1] = 0;
    }

    for (int ip = 1; ip <= grid_.npartmax; ++ip) {
        if (grid_.idp[ip - 1] <= 0) continue;

        real_t xp = grid_.xp[0 * grid_.npartmax + ip - 1];
        real_t yp = (NDIM > 1) ? grid_.xp[1 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;
        real_t zp = (NDIM > 2) ? grid_.xp[2 * grid_.npartmax + ip - 1] : 0.5 * params::boxlen;

        int level = 1;
        int icell = find_cell_by_coords(xp, yp, zp, &level);
        if (icell > 0) {
            grid_.levelp[ip - 1] = level;
            int ig = 0;
            if (icell > grid_.ncoarse) {
                ig = ((icell - grid_.ncoarse - 1) % grid_.ngridmax) + 1;
            }
            
            if (ig > 0) {
                if (grid_.numbp[ig - 1] > 0) {
                    grid_.nextp[grid_.tailp[ig - 1] - 1] = ip;
                    grid_.prevp[ip - 1] = grid_.tailp[ig - 1];
                    grid_.nextp[ip - 1] = 0;
                    grid_.tailp[ig - 1] = ip;
                    grid_.numbp[ig - 1]++;
                } else {
                    grid_.headp[ig - 1] = ip;
                    grid_.tailp[ig - 1] = ip;
                    grid_.prevp[ip - 1] = 0;
                    grid_.nextp[ip - 1] = 0;
                    grid_.numbp[ig - 1] = 1;
                }
            }
        }
    }
}

int ParticleSolver::find_cell_by_coords(real_t x, real_t y, real_t z, int* level) {
    real_t dx_coarse = params::boxlen / std::max({params::nx, params::ny, params::nz});
    int ix = std::clamp(static_cast<int>(x / dx_coarse), 0, params::nx - 1);
    int iy = std::clamp(static_cast<int>(y / dx_coarse), 0, params::ny - 1);
    int iz = std::clamp(static_cast<int>(z / dx_coarse), 0, params::nz - 1);
    int curr_cell = iz * params::nx * params::ny + iy * params::nx + ix + 1;
    
    int ilevel = 1;
    while (grid_.son[curr_cell - 1] > 0) {
        int ig = grid_.son[curr_cell - 1];
        ilevel++;
        real_t xg = grid_.xg[0 * grid_.ngridmax + ig - 1];
        real_t yg = (NDIM > 1) ? grid_.xg[1 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
        real_t zg = (NDIM > 2) ? grid_.xg[2 * grid_.ngridmax + ig - 1] : 0.5 * params::boxlen;
        int ixc = (x > xg) ? 1 : 0;
        int iyc = (y > yg) ? 1 : 0;
        int izc = (z > zg) ? 1 : 0;
        int ic = izc * 4 + iyc * 2 + ixc + 1;
        curr_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
    }
    if (level) *level = ilevel;
    return curr_cell;
}

} // namespace ramses
