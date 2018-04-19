#include "11W.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

void Clusters_Get11W() {
    // 10B with 1 extra particle. The central particle of 10B must have exactly 10 particles bonded to it.
    // The extra particle is bonded to 10B central particle but not bonded to three shell spindles of 7A in 10B

    int id_10B, spindle_10B;
    int extra_particle;

    for(id_10B=0;id_10B<n10B; id_10B++) {
        spindle_10B = hc10B[id_10B][9];
        if (num_bonds[spindle_10B] == 10) {   // s_com has 10 bonds in total (all forming the shell)

            extra_particle = get_11W_extra_particle(id_10B, spindle_10B);

            // extra particle must not be bonded to three 7A spindles in shell of 10B
            if (is_particle_bonded_to_7As(id_10B, extra_particle)) continue;

            resize_hc11W();
            populate_hc11W(id_10B, extra_particle);
            populate_s11W();
            n11W++;
        }
    }
}

void populate_hc11W(int id_10B, int extra_particle) {
    int i;
    for (i = 0; i < 10; i++){
        hc11W[n11W][i] = hc10B[id_10B][i];
    }
    hc11W[n11W][10] = extra_particle;
}

void resize_hc11W() {
    int clusSize=11;

    if (n11W == m11W) {
        hc11W = resize_2D_int(hc11W, m11W, m11W + incrStatic, clusSize, -1);
        m11W = m11W + incrStatic;
    }
}

int is_particle_bonded_to_7As(int id_10B, int extra_particle) {
    int i;

    for(i=6; i < 9; i++) {
        if (Bonds_BondCheck(extra_particle, hc10B[id_10B][i]) == 1){
            return 1;
        }
    }
    return 0;
}

int get_11W_extra_particle(int id_10B, int spindle_10B) {
    int i;
    for (i = 0; i < num_bonds[spindle_10B]; ++i) {
        if (is_particle_in_10B(bNums[spindle_10B][i], id_10B) == 0) {
            return bNums[spindle_10B][i];
        }
    }
    Error("11W extra particle not found");
    return 0;
}

int is_particle_in_10B(int particle_id, int id_10B) {
    // Returns 1 if particle_id is in id_10B else returns 0
    int i;

    for (i=0; i<9; i++) {
        if (particle_id == hc10B[id_10B][i]) {
            return 1;
        }
    }
    return 0;
}

void populate_s11W() {
    int i;

    for(i=0; i<10; i++) {
        if (s11W[hc11W[n11W][i]] == 'C') s11W[hc11W[n11W][i]] = 'B';
    }
    s11W[hc11W[n11W][10]] = 'O';
}