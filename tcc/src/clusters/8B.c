#include "8B.h"
#include "globals.h"
#include "bonds.h"
#include "tools.h"

int is_particle_in_7A(int first_7A_id, int new_particle_id);

void Clusters_Get8B() { // Detect 8B Cs clusters
    int first_7A_id;


    for (first_7A_id = 0; first_7A_id < nsp5c; ++first_7A_id) {    // loop over all 7A_i
        Clusters_8B_loop(first_7A_id, hcsp5c[first_7A_id][5], hcsp5c[first_7A_id][6]);
        Clusters_8B_loop(first_7A_id, hcsp5c[first_7A_id][6], hcsp5c[first_7A_id][5]);
    }
}

void Clusters_8B_loop(int first_7A_id, int primary_spindle, int secondary_spindle) {

    int new_particle_pointer, k, l, m;
    int new_particle_id, nbs, unc[3];
    int break_out;
    int clusSize = 8;

    for (new_particle_pointer = 0; new_particle_pointer < num_bonds[primary_spindle]; ++new_particle_pointer) {
        new_particle_id = bNums[primary_spindle][new_particle_pointer];

        if(is_particle_in_7A(first_7A_id, new_particle_id) == 0) {

            if (new_particle_id == secondary_spindle) continue; // now is new_particle_id bonded to sp5

            nbs = 0; // number of bonds
            for (k = 0; k < 5; ++k) {
                if (Bonds_BondCheck(new_particle_id, hcsp5c[first_7A_id][k])) {
                    ++nbs;
                }
            }
            if (nbs != 2) continue;

            // Now we have found the 8B Cs cluster
            if (n8B == m8B) {
                hc8B = resize_2D_int(hc8B, m8B, m8B + incrStatic, clusSize, -1);
                m8B = m8B + incrStatic;
            }

            l = 0;
            m = 3;
            break_out = 0;
            for (k = 0; k < 5; ++k) {
                if (Bonds_BondCheck(new_particle_id, hcsp5c[first_7A_id][k])) {
                    if (m == 5) {
                        break_out = 1;
                        break;
                    }
                    hc8B[n8B][m] = hcsp5c[first_7A_id][k];
                    m++;
                } else {
                    if (l == 3) {
                        break_out = 1;
                        break;
                    }
                    unc[l] = hcsp5c[first_7A_id][k];
                    l++;
                }
            }
            if (break_out == 1 || m < 5 || l < 3) continue;

            quickSort(&hc8B[n8B][3], 2);
            for (k = 0; k < 3; k++) hc8B[n8B][k] = unc[k];
            quickSort(&hc8B[n8B][0], 3);

            hc8B[n8B][5] = secondary_spindle;
            hc8B[n8B][6] = primary_spindle;
            hc8B[n8B][7] = new_particle_id;


            Cluster_Write_8B();
        }
    }
}

int is_particle_in_7A(int first_7A_id, int new_particle_id) {
    int is_in_7A = 0;
    for (int first_7A_pointer = 0; first_7A_pointer < 5; ++first_7A_pointer) {
        if (new_particle_id == hcsp5c[first_7A_id][first_7A_pointer]) {
            is_in_7A = 1;
            break;
        }
    }
    return is_in_7A;
}

void Cluster_Write_8B() {
    // hc8B key: (SP5_to_4, SP5_to_0/2, SP5_to_3, SP5_to_n1(lower), SP5_to_n1(greater), s, s_to_n1, n1)
    if (s8B[hc8B[n8B][7]] == 'C') s8B[hc8B[n8B][7]] = 'B';
    if (s8B[hc8B[n8B][0]] == 'C') s8B[hc8B[n8B][0]] = 'B';
    if (s8B[hc8B[n8B][1]] == 'C') s8B[hc8B[n8B][1]] = 'B';
    if (s8B[hc8B[n8B][2]] == 'C') s8B[hc8B[n8B][2]] = 'B';
    if (s8B[hc8B[n8B][3]] == 'C') s8B[hc8B[n8B][3]] = 'B';
    if (s8B[hc8B[n8B][4]] == 'C') s8B[hc8B[n8B][4]] = 'B';
    s8B[hc8B[n8B][5]] = 'O';
    s8B[hc8B[n8B][6]] = 'O';
    ++n8B;
}