#include "12A.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

int Clusters_Get12A() {
    // A 12A is an 11C with an extra particle bonded to only 2 other specific outer shell particles in the 11C.

    int id_11C;
    int ep;

    for(id_11C=0; id_11C<n11C; id_11C++) {
        if (num_bonds[hc11C[id_11C][0]] == 11) {

            ep = get_12A_extra_particle(id_11C);

            // Extra particle must be bonded to a specific two particles in the 11C
            if (Bonds_BondCheck(ep, hc11C[id_11C][9]) == 0 || Bonds_BondCheck(ep, hc11C[id_11C][10]) == 0) continue;

            // The extra particle should not be bonded to particles 2-8 of the 11C
            if (bond_check_12A_extra_particle(id_11C, ep) == 1) continue;

            resize_hc12A();
            populate_hc12A(id_11C, ep);
            populate_s12A();
            ++n12A;
        }
    }
}

int get_12A_extra_particle(int id_11C) {
    int i;
    // Returns id of extra particle
    // The extra particle is the one bonded to the 11C center that is not in the 11C,
    for (i = 0; i < num_bonds[hc11C[id_11C][0]]; ++i) {
        if (is_particle_in_11C(bNums[hc11C[id_11C][0]][i], id_11C) == 0) {
            return bNums[hc11C[id_11C][0]][i]; // The extra particle
        }
    }
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

int is_particle_in_11C(int particle_id, int id_11C) {
    // Return 1 if particle is in 11C, else returns 0
    int i;

    for (i=1; i<11; i++) {
        if (particle_id == hc11C[id_11C][i]) {
            return 1;
        }
    }
    return 0;
}

void populate_hc12A(int id_11C, int ep) {
    int i;

    for (i = 0; i<11; i++) {
        hc12A[n12A][i] = hc11C[id_11C][i];
    }
    hc12A[n12A][11] = ep;
}

void resize_hc12A() {
    int clusSize=12;
    if (n12A == m12A) {
        hc12A = resize_2D_int(hc12A, m12A, m12A + incrStatic, clusSize, -1);
        m12A = m12A + incrStatic;
    }
}

void populate_s12A() {
    int i;
    // hc12A key: (as 11C, extra_s)
    s12A[hc12A[n12A][0]] = 'S';
    if(s12A[hc12A[n12A][1]] != 'S') s12A[hc12A[n12A][1]] = 'O';
    if(s12A[hc12A[n12A][2]] != 'S') s12A[hc12A[n12A][2]] = 'O';
    for(i=3; i<12; i++) {
        if (s12A[hc12A[n12A][i]] == 'C') s12A[hc12A[n12A][i]] = 'B';
    }
}