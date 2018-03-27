#include "11A.h"
#include "globals.h"
#include "tools.h"
#include "bonds.h"

void Clusters_Get11A() {

    int first_6A_id, second_6A_id;
    int *first6A, *second6A;
    int scom, sother[2];

    scom = sother[0] = sother[1] = -1;

    for(first_6A_id = 0; first_6A_id < nsp4c - 1; first_6A_id++){
        first6A = hcsp4c[first_6A_id];
        for(second_6A_id = first_6A_id + 1; second_6A_id < nsp4c; second_6A_id++) {
            second6A = hcsp4c[second_6A_id];

            if(count_common_spindle_particles(first6A, second6A, &scom) == 1) {

                get_non_common_spindles(first6A, second6A, scom, sother);

                if (Check_unique_6A_rings(first6A, second6A) == 0) {
                    if (Check_6A_rings_bonded(first6A, second6A) == 1) {
                        Cluster_Write_11A(first6A, second6A, sother, scom);
                    }
                }
            }
        }
    }
}

void get_non_common_spindles(const int *first6A, const int *second6A, int scom, int *sother) {
    if (scom == first6A[4]) {
        sother[0] = first6A[5];
    }
    else {
        sother[0] = first6A[4];
    }
    if (scom == second6A[4]) {
        sother[1] = second6A[5];
    }
    else {
        sother[1] = second6A[4];
    }
}

int count_common_spindle_particles(const int *first6A, const int *second6A, int *scom) {
    int first_spindle_pointer, second_spindle_pointer, num_common_spindles;

    num_common_spindles = 0;

    for (first_spindle_pointer = 4; first_spindle_pointer < 6; first_spindle_pointer++) {
        for (second_spindle_pointer = 4; second_spindle_pointer < 6; second_spindle_pointer++) {
            if (first6A[first_spindle_pointer] == second6A[second_spindle_pointer]) {
                (*scom) = first6A[first_spindle_pointer];
                num_common_spindles++;
            }
        }
    }
    return num_common_spindles;
}

void Cluster_Write_11A(const int *first_6A, const int *second_6A, const int *sother, const int scom) {
    int clusSize=11;
    int i;

    if (n11A == m11A) {
        hc11A = resize_2D_int(hc11A, m11A, m11A + incrStatic, clusSize, -1);
        m11A = m11A + incrStatic;
    }

    hc11A[n11A][0] = first_6A[0];
    hc11A[n11A][1] = first_6A[1];
    hc11A[n11A][2] = first_6A[2];
    hc11A[n11A][3] = first_6A[3];
    hc11A[n11A][4] = second_6A[0];
    hc11A[n11A][5] = second_6A[1];
    hc11A[n11A][6] = second_6A[2];
    hc11A[n11A][7] = second_6A[3];
    hc11A[n11A][8] = sother[0];
    hc11A[n11A][9] = sother[1];
    hc11A[n11A][10] = scom;

    for(i=0; i<8; i++) {
        if (s11A[hc11A[n11A][i]] == 'C') s11A[hc11A[n11A][i]] = 'B';
    }
    for(i=8; i<10; i++) {
        if (s11A[hc11A[n11A][i]] != 'S') s11A[hc11A[n11A][i]] = 'O';
    }
    s11A[hc11A[n11A][10]] = 'S';
    ++n11A;
}

int Check_6A_rings_bonded(const int *first_6A, const int *second_6A) {
    int first_ring_pointer, second_ring_pointer, num_bonds;
    // Check if there are two bonds between each particle in ring 1 and particles in ring 2
    // Returns 1 if all ring 1 particles have 2 bonds to ring 2 particles, return 0 if not

    for(first_ring_pointer = 0; first_ring_pointer < 4; first_ring_pointer++) {
        num_bonds = 0;
        for(second_ring_pointer = 0; second_ring_pointer < 4; second_ring_pointer++) {
            if(Bonds_BondCheck(first_6A[first_ring_pointer], second_6A[second_ring_pointer])) {
                num_bonds++;
            }
        }
        if(num_bonds!=2) {
            return 0;
        }
    }
    return 1;
}

int Check_unique_6A_rings(const int *first_6A, const int *second_6A) {
    // Check that there are no common ring particles between two 6As.
    // Return 1 if the two 6A rings share a particle, else return 0
    int first_ring, second_ring;

    for(first_ring=0; first_ring<4; ++first_ring) {
        for(second_ring=0; second_ring<4; ++second_ring) {
            if(first_6A[first_ring] == second_6A[second_ring]) {
                return 1;
            }
        }
    }
    return 0;
}