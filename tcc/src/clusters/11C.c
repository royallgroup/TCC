#import "globals.h"
#import "bonds.h"
#import "tools.h"
#import "11C.h"

void Clusters_Get11C() {
    int ar[2],uncommon_spindle[2];
    int id_first_7A, id_second7A, k, l, m, ncom, common_spindle;
    int break_out;

    int first_spindle_id, first_spindle_pointer, second_7A_pointer;

    ar[0]=ar[1]=uncommon_spindle[0]=uncommon_spindle[1]=common_spindle=-1;

    for (id_first_7A=0; id_first_7A<nsp5c-1; ++id_first_7A) {
        for(first_spindle_pointer=5; first_spindle_pointer<7; first_spindle_pointer++) {
            first_spindle_id = hcsp5c[id_first_7A][first_spindle_pointer];
            for (second_7A_pointer = 0; second_7A_pointer < nmem_sp5c[first_spindle_id]; second_7A_pointer++) {
                id_second7A = mem_sp5c[first_spindle_id][second_7A_pointer];
                if(id_second7A<=id_first_7A) continue; // Dont detect the same cluster twice!

                if (get_11C_spindle_particles(uncommon_spindle, id_first_7A, id_second7A, &common_spindle) == 0) continue;

                ncom = 0;
                // need two common particles from SP5 rings

                ncom = get_bonded_7A_ring_particles(ar, id_first_7A, id_second7A, ncom);

                if (ncom != 2) continue;
                if (Bonds_BondCheck(ar[0], ar[1]) != 1) continue;

                // two common SP5 ring particles are bonded

                ncom = 0;
                for (k = 0; k < 5; ++k) {
                    if (hcsp5c[id_first_7A][k] == ar[0] || hcsp5c[id_first_7A][k] == ar[1]) continue;
                    for (l = 0; l < 5; ++l) {
                        if (hcsp5c[id_second7A][l] == ar[0] || hcsp5c[id_second7A][l] == ar[1]) continue;
                        if (Bonds_BondCheck(hcsp5c[id_first_7A][k], hcsp5c[id_second7A][l])) ++ncom;
                    }
                }
                if (ncom != 2) continue;

                // two bonds between non-common SP5 ring particles
                resize_hc11C();

                // hc11C key: (s_com, s_i, s_j, r_ca, r_cb, d_i, d_i, d_j, d_j, unc_i, unc_j)
                hc11C[n11C][0] = common_spindle;
                hc11C[n11C][1] = uncommon_spindle[0];
                hc11C[n11C][2] = uncommon_spindle[1];
                hc11C[n11C][3] = ar[0];
                hc11C[n11C][4] = ar[1];

                l = 5;
                m = 7;
                break_out = 0;
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], ar[0]) && hcsp5c[id_first_7A][k] != ar[1]) {
                        if (l == 7) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][l] = hcsp5c[id_first_7A][k];
                        l++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], ar[1]) && hcsp5c[id_first_7A][k] != ar[0]) {
                        if (l == 7) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][l] = hcsp5c[id_first_7A][k];
                        l++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], ar[0]) && hcsp5c[id_second7A][k] != ar[1]) {
                        if (m == 9) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][m] = hcsp5c[id_second7A][k];
                        m++;
                    }
                }
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], ar[1]) && hcsp5c[id_second7A][k] != ar[0]) {
                        if (m == 9) {
                            break_out = 1;
                            break;
                        }
                        hc11C[n11C][m] = hcsp5c[id_second7A][k];
                        m++;
                    }
                }
                if (break_out == 1 || l < 7 || m < 9) continue;

                // Check that the bonded non-common particles are bonded
                if (Bonds_BondCheck(hc11C[n11C][5], hc11C[n11C][7]) == 0) continue;
                if (Bonds_BondCheck(hc11C[n11C][6], hc11C[n11C][8]) == 0) continue;

                // Get the ID's of the non-common particles
                for (k = 0; k < 5; ++k) {
                    if (Bonds_BondCheck(hcsp5c[id_first_7A][k], hc11C[n11C][5]) &&
                        Bonds_BondCheck(hcsp5c[id_first_7A][k], hc11C[n11C][6])) {
                        hc11C[n11C][9] = hcsp5c[id_first_7A][k];
                    }
                    if (Bonds_BondCheck(hcsp5c[id_second7A][k], hc11C[n11C][7]) &&
                        Bonds_BondCheck(hcsp5c[id_second7A][k], hc11C[n11C][8])) {
                        hc11C[n11C][10] = hcsp5c[id_second7A][k];
                    }
                }
                quickSort(&hc11C[n11C][1], 2);
                quickSort(&hc11C[n11C][3], 2);
                quickSort(&hc11C[n11C][5], 4);
                quickSort(&hc11C[n11C][9], 2);

                Cluster_Write_11C();

                ++n11C;
            }
        }
    }
}

int get_11C_spindle_particles(int *uncommon_spindle, int id_first_7A, int id_second7A, int *common_spindle) {
    int num_common_spindles = 0;

    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][5] == hcsp5c[id_second7A][6]) {
        (*common_spindle) = hcsp5c[id_first_7A][5];
        uncommon_spindle[0] = hcsp5c[id_first_7A][6];
        uncommon_spindle[1] = hcsp5c[id_second7A][5];
        ++num_common_spindles;
    }
    if (hcsp5c[id_first_7A][6] == hcsp5c[id_second7A][5]) {
        (*common_spindle) = hcsp5c[id_first_7A][6];
        uncommon_spindle[0] = hcsp5c[id_first_7A][5];
        uncommon_spindle[1] = hcsp5c[id_second7A][6];
        ++num_common_spindles;
    }

    if (num_common_spindles == 1) return 1;
    else return 0;
}

int get_bonded_7A_ring_particles(int *ar, int id_first_7A, int id_second7A, int ncom) {
    int first_ring_pointer, second_ring_pointer;

    for (first_ring_pointer = 0; first_ring_pointer < 5; ++first_ring_pointer) {
        for (second_ring_pointer = 0; second_ring_pointer < 5; ++second_ring_pointer) {
            if (hcsp5c[id_first_7A][first_ring_pointer] == hcsp5c[id_second7A][second_ring_pointer]) {
                if (ncom == 2) {
                    ++ncom;
                    break;
                }
                ar[ncom++] = hcsp5c[id_first_7A][first_ring_pointer];
                break;
            }
        }
        if (ncom > 2) break;
    }
    return ncom;
}

void resize_hc11C() {
    int clusSize=11;

    if(n11C == m11C) {
        hc11C=resize_2D_int(hc11C,m11C,m11C+incrStatic,clusSize,-1);
        m11C=m11C+incrStatic;
    }
}

void Cluster_Write_11C() {
    int i;

    s11C[hc11C[n11C][0]] = 'S';
    if(s11C[hc11C[n11C][1]] != 'S') s11C[hc11C[n11C][1]] = 'O';
    if(s11C[hc11C[n11C][2]] != 'S') s11C[hc11C[n11C][2]] = 'O';
    for(i=3; i< 11; i++) {
        if (s11C[hc11C[n11C][i]] == 'C') s11C[hc11C[n11C][i]] = 'B';
    }
}