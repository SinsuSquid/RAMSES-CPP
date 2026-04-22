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

void LoadBalancer::compute_new_cpu_map(std::vector<int>& cpu_map_new) {
    auto& mpi = MpiManager::instance();
    int cells_per_cpu = grid_.ncell / mpi.size();
    for (int i = 1; i <= grid_.ncell; ++i) {
        int target_rank = (i - 1) / cells_per_cpu;
        cpu_map_new[i] = std::min(target_rank + 1, mpi.size());
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
                    for (int iv = 1; iv <= 5; ++iv) p.uold[ic - 1][iv - 1] = grid_.uold(ind_cell, iv);
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
                for (int iv = 1; iv <= 5; ++iv) grid_.uold(ind_cell, iv) = p.uold[ic-1][iv-1];
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
