#include "simple_cluster_methods.h"
#include "12K.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

//!  A 12K clusters is an 11A with an extra particle bonded to three of the ring particles of the 11A
/*!
*  Find 12K clusters
*  A 12K is constructed of from an 11A and an extra particle where:
*      - The additional particle is bonded to three mutually bonded sp4 ring particles of the 11A cluster.
*
*      Since there are 8 possible sites for bonding on each 11A there may be multiple 12K for each 11A
*
*  Cluster output: BBBBBBBBOOSB
*  Storage order: as_for_11A x 11, extra_particle
*
*/
void Clusters_Get12K() {
    int pointer_11A, ring_number;
    int sp3_rings[8][3] = {-1};
    int *sp3_ring;

    // Loop over all 11A clusters
    for(pointer_11A = 0; pointer_11A < n11A; pointer_11A++) {
        int *cluster_11A = hc11A[pointer_11A];
        // Find the particle IDs of the 8 sp3 rings made from the rings of 11A
        get_12K_ring_bonds(cluster_11A, sp3_rings);

        // Loop through the 8 sp3 rings made by the sp4s of 11A and see if a particle is attached to any of them
        for(ring_number = 0; ring_number < 8; ring_number++) {
            sp3_ring = sp3_rings[ring_number];
            find_12K_cluster(cluster_11A, sp3_ring);
        }
    }
}

void find_12K_cluster(int *parent_11A_cluster, const int *sp3_ring) {
    // Check if there is a particle attached to the sp3 ring of an 11A, if there is write a 12K
    // It is possible that there are two particles attached to the ring particles in which case we ignore this cluster
    int ep = -1;

    int num_attached_particles = 0;
    // loop through all particles bonded to the first sp3 ring particle
    for (int i = 0; i < num_bonds[sp3_ring[0]]; i++) {
        int particle_id = bond_list[sp3_ring[0]][i];
        if (Bonds_BondCheck(sp3_ring[1], particle_id) == 1) {
            if (Bonds_BondCheck(sp3_ring[2], particle_id) == 1) {
                if (is_particle_in_cluster(parent_11A_cluster, 11, particle_id) == 0) {
                    num_attached_particles++;
                    ep = particle_id;
                }
            }
        }
    }
    if (num_attached_particles == 1) {
        Cluster_Write_12K(ep, parent_11A_cluster);
    }
}

void get_12K_ring_bonds(int *cluster_11A, int (*sp3_rings)[3]) {

    for(int first_sp4_ring_pointer = 0; first_sp4_ring_pointer < 4; first_sp4_ring_pointer++) {
        int m = 0;
        sp3_rings[first_sp4_ring_pointer][m] = cluster_11A[first_sp4_ring_pointer];
        for(int second_sp4_ring_pointer = 4; second_sp4_ring_pointer < 8; second_sp4_ring_pointer++) {
            if (Bonds_BondCheck(cluster_11A[first_sp4_ring_pointer], cluster_11A[second_sp4_ring_pointer])) {
                m++;
                sp3_rings[first_sp4_ring_pointer][m] = cluster_11A[second_sp4_ring_pointer];
            }
        }
    }

    for(int second_sp4_ring_pointer = 4; second_sp4_ring_pointer < 8; second_sp4_ring_pointer++) {
        int m = 0;
        sp3_rings[second_sp4_ring_pointer][m] = cluster_11A[second_sp4_ring_pointer];
        for(int first_sp4_ring_pointer = 0; first_sp4_ring_pointer < 4; first_sp4_ring_pointer++) {
            if (Bonds_BondCheck(cluster_11A[second_sp4_ring_pointer], cluster_11A[first_sp4_ring_pointer])) {
                m++;
                sp3_rings[second_sp4_ring_pointer][m] = cluster_11A[first_sp4_ring_pointer];
            }
        }
    }
}

void Cluster_Write_12K(int ep, const int *parent_11A_cluster) {
    int clusSize=12;

    if(n12K == m12K) {
        hc12K=resize_2D_int(hc12K,m12K,m12K+incrStatic,clusSize,-1);
        m12K=m12K+incrStatic;
    }

    for (int i = 0; i < 11; i++) hc12K[n12K][i] = parent_11A_cluster[i];
    hc12K[n12K][11]=ep;

    for (int i = 0; i < 8; i++) {
        if (s12K[hc12K[n12K][i]] == 'C') s12K[hc12K[n12K][i]] = 'B';
    }
    if(s12K[hc12K[n12K][8]] != 'S') s12K[hc12K[n12K][8]] = 'O';
    if(s12K[hc12K[n12K][9]] != 'S') s12K[hc12K[n12K][9]] = 'O';
    s12K[hc12K[n12K][10]] = 'S';
    if(s12K[hc12K[n12K][11]] == 'C') s12K[hc12K[n12K][11]] = 'B';
    n12K++;
}