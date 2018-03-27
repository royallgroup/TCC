#include "12K.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get12K() {
    // 12K clusters are an 11A with an extra particle bonded to three of the ring particles of the 11A
    // Since there are 8 possible places to attach an extra particle, there may be more than 1 12K for each 11A.

    int ptr_11A, ring_number;
    int sp3_rings[8][3] = {-1};
    int *sp3_ring;

    // Loop over all 11A clusters
    for(ptr_11A=0; ptr_11A<n11A; ptr_11A++) {

        // Find the particle IDs of the 8 sp3 rings made from the rings of 11A
        get_12K_ring_bonds(ptr_11A, sp3_rings);

        // Loop through the 8 sp3 rings made by the sp4s of 11A and see if a particle is attached to any of them
        for(ring_number=0; ring_number<8; ring_number++) {
            sp3_ring = sp3_rings[ring_number];
            find_12K_cluster(ptr_11A, sp3_ring);
        }
    }
}

void find_12K_cluster(int ptr_11A, const int *sp3_ring) {
    // Check if there is a particle attached to the sp3 ring of an 11A, if there is write a 12K
    // It is possible that there are two particles attached to the ring particles in which case we ignore this cluster
    int i, num_attached_particles;
    int particle_id, ep = {-1};

    num_attached_particles = 0;
    // loop through all particles bonded to the first sp3 ring particle
    for (i = 0; i < num_bonds[sp3_ring[0]]; i++) {
        particle_id = bNums[sp3_ring[0]][i];
        if (Bonds_BondCheck(sp3_ring[1], particle_id) == 1) {
            if (Bonds_BondCheck(sp3_ring[2], particle_id) == 1) {
                if (is_particle_in_11A(ptr_11A, particle_id) == 0) {
                    num_attached_particles++;
                    ep = particle_id;
                }
            }
        }
    }
    if (num_attached_particles == 1) {
        Cluster_Write_12K(ep, ptr_11A);
    }
}

int is_particle_in_11A(int id_11A, int id_particle) {
    // Returns 0 if particle is not in specified 11A
    int i;
    for (i = 0; i < 11; i++) {
        if (id_particle == hc11A[id_11A][i]) {
            return 1;
        }
    }
    return 0;
}

void get_12K_ring_bonds(int id_11A, int (*sp3_rings)[3]) {

    int m;
    int particle_1, particle_2;

    for(particle_1=0; particle_1<4; particle_1++) {
        m = 0;
        sp3_rings[particle_1][m] = hc11A[id_11A][particle_1];
        for(particle_2=4; particle_2<8; particle_2++) {
            if (Bonds_BondCheck(hc11A[id_11A][particle_1], hc11A[id_11A][particle_2])) {
                m++;
                sp3_rings[particle_1][m] = hc11A[id_11A][particle_2];
            }
        }
    }
    for(particle_1=4; particle_1<8; particle_1++) {
        m = 0;
        sp3_rings[particle_1][m] = hc11A[id_11A][particle_1];
        for(particle_2=0; particle_2<4; particle_2++) {
            if (Bonds_BondCheck(hc11A[id_11A][particle_1], hc11A[id_11A][particle_2])) {
                m++;
                sp3_rings[particle_1][m] = hc11A[id_11A][particle_2];
            }
        }
    }
}

void Cluster_Write_12K(int ep, int id_11A) {

    int i;
    int clusSize=12;

    if(n12K == m12K) {
        hc12K=resize_2D_int(hc12K,m12K,m12K+incrStatic,clusSize,-1);
        m12K=m12K+incrStatic;
    }
    // hc12K key: (SP4 going up, sd going up, scom, ep)

    for (i=0; i<11; i++) hc12K[n12K][i] = hc11A[id_11A][i];
    hc12K[n12K][11]=ep;

    for (i=0; i<8; i++) {
        if (s12K[hc12K[n12K][i]] == 'C') s12K[hc12K[n12K][i]] = 'B';
    }
    if(s12K[hc12K[n12K][8]] != 'S') s12K[hc12K[n12K][8]] = 'O';
    if(s12K[hc12K[n12K][9]] != 'S') s12K[hc12K[n12K][9]] = 'O';
    s12K[hc12K[n12K][10]] = 'S';
    if(s12K[hc12K[n12K][11]] == 'C') s12K[hc12K[n12K][11]] = 'B';
    n12K++;
}