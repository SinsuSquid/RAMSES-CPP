#include "ramses/LoadBalancer.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Hilbert.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <algorithm>
#include <map>
#include <cmath>
#include <cstring>

#ifdef RAMSES_USE_MPI
#include <mpi.h>
#endif

namespace ramses {

void LoadBalancer::balance() {
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) return;

    std::vector<int> cpu_map_new(grid_.ncell + 1);
    compute_new_cpu_map(cpu_map_new);
    move_grids(cpu_map_new);
}

void LoadBalancer::calculate_hilbert_keys() {
    int n_res = 1 << grid_.nlevelmax;
    int myid = MpiManager::instance().rank() + 1;

    auto calc_for_cell = [&](int i) {
        real_t xc[3];
        grid_.get_cell_center(i, xc);
        
        std::vector<int> ix(1), iy(1), iz(1);
        ix[0] = std::min(n_res - 1, static_cast<int>(std::floor(xc[0] * n_res)));
        if (NDIM > 1) iy[0] = std::min(n_res - 1, static_cast<int>(std::floor(xc[1] * n_res)));
        if (NDIM > 2) iz[0] = std::min(n_res - 1, static_cast<int>(std::floor(xc[2] * n_res)));
        
        std::vector<qdp_t> order;
        if (NDIM == 1) Hilbert::hilbert1d(ix, order);
        else if (NDIM == 2) Hilbert::hilbert2d(ix, iy, order, grid_.nlevelmax);
        else Hilbert::hilbert3d(ix, iy, iz, order, grid_.nlevelmax);
        
        grid_.hilbert_keys[i-1] = order[0];
    };

    // 1. Coarse cells
    for (int i = 1; i <= grid_.ncoarse; ++i) calc_for_cell(i);

    // 2. Fine cells
    for (int il = 1; il <= grid_.nlevelmax; ++il) {
        int ig = grid_.headl(myid, il);
        while (ig > 0) {
            for (int ic = 1; ic <= constants::twotondim; ++ic) {
                int id = grid_.ncoarse + (ic - 1) * grid_.ngridmax + ig;
                calc_for_cell(id);
            }
            ig = grid_.next[ig - 1];
        }
    }
}

void LoadBalancer::compute_new_cpu_map(std::vector<int>& cpu_map_new) {
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) {
        std::fill(cpu_map_new.begin(), cpu_map_new.end(), 1);
        return;
    }

    struct CellKey {
        int index;
        qdp_t key;
    };
    std::vector<CellKey> leaf_cells;
    for (int i = 1; i <= grid_.ncell; ++i) {
        if (grid_.son[i-1] == 0) leaf_cells.push_back({i, grid_.hilbert_keys[i-1]});
    }

    std::sort(leaf_cells.begin(), leaf_cells.end(), [](const CellKey& a, const CellKey& b) {
        return a.key < b.key;
    });

    int nleaf = leaf_cells.size();
    int leaves_per_cpu = nleaf / mpi.size();
    int extra_leaves = nleaf % mpi.size();

    int current_leaf = 0;
    for (int icpu = 1; icpu <= mpi.size(); ++icpu) {
        int count = leaves_per_cpu + (icpu <= extra_leaves ? 1 : 0);
        for (int i = 0; i < count; ++i) {
            cpu_map_new[leaf_cells[current_leaf].index - 1] = icpu;
            current_leaf++;
        }
    }

    for (int ilevel = grid_.nlevelmax; ilevel >= 1; --ilevel) {
        for (int icpu = 1; icpu <= mpi.size(); ++icpu) {
            int igrid = grid_.headl(icpu, ilevel);
            while (igrid > 0) {
                int father_cell = grid_.father[igrid - 1];
                if (father_cell > 0) {
                    int first_child = grid_.ncoarse + (1 - 1) * grid_.ngridmax + igrid;
                    cpu_map_new[father_cell-1] = cpu_map_new[first_child-1];
                }
                igrid = grid_.next[igrid - 1];
            }
        }
    }
}

void LoadBalancer::remove_grid_from_list(int igrid, int ilevel, int icpu) {
    int prev_idx = grid_.prev[igrid - 1];
    int next_idx = grid_.next[igrid - 1];

    if (prev_idx > 0) grid_.next[prev_idx - 1] = next_idx;
    else grid_.headl(icpu, ilevel) = next_idx;

    if (next_idx > 0) grid_.prev[next_idx - 1] = prev_idx;
    else grid_.taill(icpu, ilevel) = prev_idx;

    grid_.numbl(icpu, ilevel)--;
}

void LoadBalancer::move_grids(const std::vector<int>& cpu_map_new) {
    auto& mpi = MpiManager::instance();
    int myid = mpi.rank() + 1;

    std::map<int, std::vector<OctPacket>> send_queues;
    std::vector<int> grids_to_remove;

    int ncells_oct = (1 << NDIM);

    for (int ilevel = 1; ilevel <= grid_.nlevelmax; ++ilevel) {
        int igrid = grid_.headl(myid, ilevel);
        while (igrid > 0) {
            int next_grid = grid_.next[igrid - 1];
            int father_cell = grid_.father[igrid - 1];
            int target_cpu = cpu_map_new[father_cell-1];

            if (target_cpu != myid) {
                OctPacket p;
                std::memset(&p, 0, sizeof(OctPacket));
                p.ilevel = ilevel;
                grid_.get_cell_center(father_cell, p.father_x);

                for (int d = 0; d < 3; ++d) p.xg[d] = (d < NDIM) ? grid_.xg[d * grid_.ngridmax + igrid - 1] : 0.0;
                for (int n = 0; n < 6; ++n) p.nbor[n] = (n < 2*NDIM) ? grid_.nbor[n * grid_.ngridmax + igrid - 1] : 0;
                
                for (int ic = 1; ic <= ncells_oct; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    for (int iv = 1; iv <= std::min(grid_.nvar, 30); ++iv) {
                        p.uold[(ic - 1) * 30 + (iv - 1)] = grid_.uold(idc, iv);
                        p.unew[(ic - 1) * 30 + (iv - 1)] = grid_.unew(idc, iv);
                    }
                    p.hilbert_keys[ic - 1] = grid_.hilbert_keys[idc-1];
                }
                send_queues[target_cpu].push_back(p);
                grids_to_remove.push_back(igrid);
                remove_grid_from_list(igrid, ilevel, myid);
                grid_.son[father_cell-1] = 0; // Clear son pointer on source
            }
            igrid = next_grid;
        }
    }

#ifdef RAMSES_USE_MPI
    int ncpu = mpi.size();
    std::vector<int> send_counts(ncpu, 0), recv_counts(ncpu, 0);
    for (auto const& [cpu, queue] : send_queues) send_counts[cpu-1] = queue.size();

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<MPI_Request> requests;
    std::vector<std::vector<OctPacket>> recv_queues(ncpu);

    for (int i = 0; i < ncpu; ++i) {
        if (i == mpi.rank()) continue;
        if (send_counts[i] > 0) {
            MPI_Request req;
            MPI_Isend(send_queues[i+1].data(), send_counts[i] * sizeof(OctPacket), MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
        if (recv_counts[i] > 0) {
            recv_queues[i].resize(recv_counts[i]);
            MPI_Request req;
            MPI_Irecv(recv_queues[i].data(), recv_counts[i] * sizeof(OctPacket), MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    for (int il = 1; il <= grid_.nlevelmax; ++il) {
        for (int i = 0; i < ncpu; ++i) {
            for (const auto& p : recv_queues[i]) {
                if (p.ilevel != il) continue;
                
                int igrid = grid_.get_free_grid();
                if (igrid == 0) {
                    std::cerr << "[LoadBalancer] Error: No free grids available during migration!" << std::endl;
                    continue;
                }

                int father_cell = grid_.find_cell_by_coords(p.father_x, p.ilevel - 1);
                grid_.father[igrid - 1] = father_cell;
                grid_.son[father_cell-1] = igrid;

                for(int d=0; d<NDIM; ++d) grid_.xg[d * grid_.ngridmax + igrid - 1] = p.xg[d];
                for(int n=0; n<2*NDIM; ++n) grid_.nbor[n * grid_.ngridmax + igrid - 1] = p.nbor[n];
                
                for (int ic = 1; ic <= ncells_oct; ++ic) {
                    int idc = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    for (int iv = 1; iv <= std::min(grid_.nvar, 30); ++iv) {
                        grid_.uold(idc, iv) = p.uold[(ic - 1) * 30 + (iv - 1)];
                        grid_.unew(idc, iv) = p.unew[(ic - 1) * 30 + (iv - 1)];
                    }
                    grid_.hilbert_keys[idc-1] = p.hilbert_keys[ic - 1];
                    grid_.cpu_map[idc-1] = myid;
                }

                grid_.add_to_level_list(igrid, p.ilevel);
            }
        }
    }
#endif

    for (int i = 1; i <= grid_.ncell; ++i) grid_.cpu_map[i-1] = cpu_map_new[i-1];

    for (int igrid : grids_to_remove) {
        grid_.free_grid(igrid);
    }
}

} // namespace ramses
