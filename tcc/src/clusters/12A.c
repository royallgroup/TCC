#include <clusters/simple_cluster_methods.h>
#include "12A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get12A() {
    // A 12A is an 11C with an extra particle bonded to only 2 other specific outer shell particles in the 11C.

    int first_11C_id;
    int extra_particle;

    for(first_11C_id = 0; first_11C_id < n11C; first_11C_id++) {
        int *first_11C_cluster = hc11C[first_11C_id];
        if (num_bonds[first_11C_cluster[0]] == 11) {

            extra_particle = get_12A_extra_particle(first_11C_cluster);

            // Extra particle must be bonded to a specific two particles in the 11C
            if (Bonds_BondCheck(extra_particle, first_11C_cluster[9]) == 1 && Bonds_BondCheck(extra_particle, first_11C_cluster[10]) == 1) {

                // The extra particle should not be bonded to particles 2-8 of the 11C
                if (bond_check_12A_extra_particle(first_11C_cluster, extra_particle) == 0) {

                    Write_12A(first_11C_cluster, extra_particle);
                }
            }
        }
    }
}

int get_12A_extra_particle(int *parent_11C_cluster) {
    // Returns id of extra particle
    // The extra particle is the one bonded to the 11C center that is not in the 11C,
    for (int i = 0; i < num_bonds[parent_11C_cluster[0]]; ++i) {
        int extra_particle = bNums[parent_11C_cluster[0]][i];
        if (is_particle_in_cluster(parent_11C_cluster, 11, extra_particle) == 0) {
            return extra_particle;
        }
    }
    Error("12A extra particle not found.");
    return 0;
}

int bond_check_12A_extra_particle(int *first_11C_cluster, int extra_particle) {
    // Return 1 if particle is bonded to particles 1-8 of the 11C, else return 0
    for (int i = 1; i < 9; ++i) {
        if (Bonds_BondCheck(extra_particle, first_11C_cluster[i])){
            return 1;
        }
    }
    return 0;
}

void Write_12A(const int *first_11C_cluster, int ep) {
    int clusSize=12;

    if (n12A == m12A) {
        hc12A = resize_2D_int(hc12A, m12A, m12A + incrStatic, clusSize, -1);
        m12A = m12A + incrStatic;
    }

    for (int i = 0; i < 11; i++) {
        hc12A[n12A][i] = first_11C_cluster[i];
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
