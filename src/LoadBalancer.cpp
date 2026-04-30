#include "ramses/LoadBalancer.hpp"
#include "ramses/MpiManager.hpp"
#include <iostream>
#include <algorithm>
#include <map>

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
    for (int i = 1; i <= grid_.ncell; ++i) {
        real_t xc[3];
        grid_.get_cell_center(i, xc);
        
        std::vector<int> ix(1), iy(1), iz(1);
        ix[0] = static_cast<int>(std::floor(xc[0] * n_res));
        if (NDIM > 1) iy[0] = static_cast<int>(std::floor(xc[1] * n_res));
        if (NDIM > 2) iz[0] = static_cast<int>(std::floor(xc[2] * n_res));
        
        std::vector<qdp_t> order;
        if (NDIM == 1) Hilbert::hilbert1d(ix, order);
        else if (NDIM == 2) Hilbert::hilbert2d(ix, iy, order, grid_.nlevelmax);
        else Hilbert::hilbert3d(ix, iy, iz, order, grid_.nlevelmax);
        
        grid_.hilbert_keys[i] = order[0];
    }
}

void LoadBalancer::compute_new_cpu_map(std::vector<int>& cpu_map_new) {
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) {
        std::fill(cpu_map_new.begin(), cpu_map_new.end(), 1);
        return;
    }

    // 1. Collect all leaf cells
    struct CellKey {
        int index;
        qdp_t key;
    };
    std::vector<CellKey> leaf_cells;
    for (int i = 1; i <= grid_.ncell; ++i) {
        if (grid_.son[i] == 0) leaf_cells.push_back({i, grid_.hilbert_keys[i]});
    }

    // 2. Sort by Hilbert key
    std::sort(leaf_cells.begin(), leaf_cells.end(), [](const CellKey& a, const CellKey& b) {
        return a.key < b.key;
    });

    // 3. Partition into ncpu segments
    int nleaf = leaf_cells.size();
    int leaves_per_cpu = nleaf / mpi.size();
    int extra_leaves = nleaf % mpi.size();

    int current_leaf = 0;
    for (int icpu = 1; icpu <= mpi.size(); ++icpu) {
        int count = leaves_per_cpu + (icpu <= extra_leaves ? 1 : 0);
        for (int i = 0; i < count; ++i) {
            cpu_map_new[leaf_cells[current_leaf].index] = icpu;
            current_leaf++;
        }
    }

    // 4. Propagate assignments up the tree (coarse cells)
    for (int ilevel = grid_.nlevelmax; ilevel >= 1; --ilevel) {
        for (int icpu = 1; icpu <= mpi.size(); ++icpu) {
            int igrid = grid_.headl(icpu, ilevel);
            while (igrid > 0) {
                int father_cell = grid_.father[igrid - 1];
                if (father_cell > 0) {
                    // Assign father to the CPU of its first child
                    int first_child = grid_.ncoarse + (1 - 1) * grid_.ngridmax + igrid;
                    cpu_map_new[father_cell] = cpu_map_new[first_child];
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

    // 1. Pack outgoing grids
    for (int ilevel = 1; ilevel <= grid_.nlevelmax; ++ilevel) {
        int igrid = grid_.headl(myid, ilevel);
        while (igrid > 0) {
            int next_grid = grid_.next[igrid - 1];
            int father_cell = grid_.father[igrid - 1];
            int target_cpu = cpu_map_new[father_cell];

            if (target_cpu != myid) {
                OctPacket p;
                p.ilevel = ilevel;
                p.father = father_cell;
                for (int d = 0; d < 3; ++d) p.xg[d] = grid_.get_xg(igrid, d + 1);
                for (int n = 0; n < 6; ++n) p.nbor[n] = grid_.get_nbor(igrid, n + 1);
                for (int ic = 1; ic <= 8; ++ic) {
                    int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                    for (int iv = 1; iv <= grid_.nvar; ++iv) {
                        p.uold[ic - 1][iv - 1] = grid_.uold(ind_cell, iv);
                        p.unew[ic - 1][iv - 1] = grid_.unew(ind_cell, iv);
                    }
                    p.hilbert_keys[ic - 1] = grid_.hilbert_keys[ind_cell];
                }
                send_queues[target_cpu].push_back(p);
                grids_to_remove.push_back(igrid);
                remove_grid_from_list(igrid, ilevel, myid);
            }
            igrid = next_grid;
        }
    }

    // 2. Perform MPI Exchange
#ifdef RAMSES_USE_MPI
    int ncpu = mpi.size();
    std::vector<int> send_counts(ncpu, 0), recv_counts(ncpu, 0);
    for (auto const& [rank, queue] : send_queues) send_counts[rank-1] = queue.size();

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

    // 3. Unpack incoming grids and reconstruct tree
    for (int i = 0; i < ncpu; ++i) {
        for (const auto& p : recv_queues[i]) {
            if (grid_.numbf <= 0) continue;
            int igrid = grid_.headf;
            grid_.headf = grid_.next[igrid - 1];
            if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = 0;
            grid_.numbf--;

            // Reconstruct connection
            grid_.father[igrid - 1] = p.father;
            grid_.son[p.father] = igrid;
            for(int d=0; d<3; ++d) grid_.get_xg(igrid, d+1) = p.xg[d];
            for(int n=0; n<6; ++n) grid_.get_nbor(igrid, n+1) = p.nbor[n];
            
            // Unpack Hydro
            for (int ic = 1; ic <= 8; ++ic) {
                int ind_cell = grid_.ncoarse + (ic - 1) * grid_.ngridmax + igrid;
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    grid_.uold(ind_cell, iv) = p.uold[ic-1][iv-1];
                    grid_.unew(ind_cell, iv) = p.unew[ic-1][iv-1];
                }
                grid_.hilbert_keys[ind_cell] = p.hilbert_keys[ic - 1];
                grid_.cpu_map[ind_cell] = myid;
            }

            // Attach to local list
            int head = grid_.headl(myid, p.ilevel);
            grid_.next[igrid - 1] = head;
            if (head > 0) grid_.prev[head - 1] = igrid;
            grid_.headl(myid, p.ilevel) = igrid;
            grid_.prev[igrid - 1] = 0;
            grid_.numbl(myid, p.ilevel)++;
        }
    }
#endif

    // 3. Update local maps
    for (int i = 1; i <= grid_.ncell; ++i) grid_.cpu_map[i] = cpu_map_new[i];

    // 4. Return removed grids to free list
    for (int igrid : grids_to_remove) {
        grid_.next[igrid - 1] = grid_.headf;
        if (grid_.headf > 0) grid_.prev[grid_.headf - 1] = igrid;
        grid_.headf = igrid;
        grid_.prev[igrid - 1] = 0;
        grid_.numbf++;
    }
}

} // namespace ramses
