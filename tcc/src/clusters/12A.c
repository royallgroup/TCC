#include <clusters/simple_cluster_methods.h>
#include "12A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get12A() {
    // A 12A is an 11C with an extra particle bonded to only 2 other specific outer shell particles in the 11C.

    int parent_11C_id;
    int ep;

    for(parent_11C_id = 0; parent_11C_id < n11C; parent_11C_id++) {
        int *parent_11C_cluster = hc11C[parent_11C_id];
        if (num_bonds[parent_11C_cluster[0]] == 11) {

            ep = get_12A_extra_particle(parent_11C_cluster);

            // Extra particle must be bonded to a specific two particles in the 11C
            if (Bonds_BondCheck(ep, parent_11C_cluster[9]) == 0 || Bonds_BondCheck(ep, parent_11C_cluster[10]) == 0) continue;

            // The extra particle should not be bonded to particles 2-8 of the 11C
            if (bond_check_12A_extra_particle(parent_11C_id, ep) == 1) continue;

            Write_12A(parent_11C_id, ep);
        }
    }
}

int get_12A_extra_particle(int *parent_11C_cluster) {
    int i;
    // Returns id of extra particle
    // The extra particle is the one bonded to the 11C center that is not in the 11C,
    for (i = 0; i < num_bonds[parent_11C_cluster[0]]; ++i) {
        int extra_particle = bNums[parent_11C_cluster[0]][i];
        if (is_particle_in_cluster(parent_11C_cluster, 11, extra_particle) == 0) {
            return extra_particle; // The extra particle
        }
    }
    Error("12A extra particle not found.");
    return 0;
}

int bond_check_12A_extra_particle(int id_11C, int extra_particle) {
    // Return 1 if particle is bonded to particles 1-8 of the 11C, else return 0
    int i;

    for (i = 1; i < 9; ++i) {
        if (Bonds_BondCheck(extra_particle, hc11C[id_11C][i])){
            return 1;
        }
    }
    return 0;
}

void Write_12A(int id_11C, int ep) {
    int clusSize=12;

    if (n12A == m12A) {
        hc12A = resize_2D_int(hc12A, m12A, m12A + incrStatic, clusSize, -1);
        m12A = m12A + incrStatic;
    }

    for (int i = 0; i < 11; i++) {
        hc12A[n12A][i] = hc11C[id_11C][i];
    }
    hc12A[n12A][11] = ep;

    // hc12A key: (as 11C, extra_s)
    s12A[hc12A[n12A][0]] = 'S';
    if(s12A[hc12A[n12A][1]] != 'S') s12A[hc12A[n12A][1]] = 'O';
    if(s12A[hc12A[n12A][2]] != 'S') s12A[hc12A[n12A][2]] = 'O';
    for(int i = 3; i < 12; i++) {
        if (s12A[hc12A[n12A][i]] == 'C') s12A[hc12A[n12A][i]] = 'B';
    }

    n12A++;
}
