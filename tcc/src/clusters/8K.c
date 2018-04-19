#include "simple_cluster_methods.h"
#include "8K.h"
#include "globals.h"
#include "tools.h"

int is_particle_in_5A_ring(const int *second_5A_cluster, int first_5A_ring_particle, int *cp);

void Clusters_Get8K() {
    int k, l, m, n;
    int cp[2], unc[3], scom, sother[2];
    int clusSize = 8;

    cp[0] = cp[1] = unc[0] = unc[1] = unc[2] = scom = sother[0] = sother[1] = -1;

    for (int first_5A_id = 0; first_5A_id < nsp3c; ++first_5A_id) {
        int *first_5A_cluster = hcsp3c[first_5A_id];
        for (int first_5A_ring_pointer = 0; first_5A_ring_pointer < 3; first_5A_ring_pointer++) {
            int first_5A_ring_particle = first_5A_cluster[first_5A_ring_pointer];
            for (int second_5A_pointer = 0;
                 second_5A_pointer < nmem_sp3c[first_5A_ring_particle]; ++second_5A_pointer) {
                int second_5A_id = mem_sp3c[first_5A_ring_particle][second_5A_pointer];
                if (second_5A_id <= first_5A_id) continue;
                int *second_5A_cluster = hcsp3c[second_5A_id];

                if (is_particle_in_5A_ring(second_5A_cluster, first_5A_ring_particle, cp) == 1) {

                    if (count_common_spindles_between_5As(first_5A_cluster, second_5A_cluster, &scom) == 1) {

                        m = 0;
                        for (k = 0; k < 3; ++k) {
                            if (k != first_5A_ring_pointer) {
                                if (first_5A_cluster[k] >= first_5A_ring_particle) {
                                    for (l = 0; l < 3; ++l) {
                                        if (first_5A_cluster[k] == second_5A_cluster[l]) {
                                            cp[1] = first_5A_cluster[k];
                                            m++;
                                        }
                                    }
                                }
                            }
                        }
                        if (m != 1) continue;

                        if (first_5A_cluster[3] == scom) {
                            sother[0] = first_5A_cluster[4];
                        } else {
                            sother[0] = first_5A_cluster[3];
                        }
                        if (second_5A_cluster[3] == scom) {
                            sother[1] = second_5A_cluster[4];
                        } else {
                            sother[1] = second_5A_cluster[3];
                        }

                        k = 0;    // find particle from SP3 ring of sp3c_i which is common with 5A mem_sp3c[sp3c[i][j2]][j]
                        for (l = 0; l < 3; ++l) {
                            if (first_5A_cluster[l] == cp[0] || first_5A_cluster[l] == cp[1]) continue;
                            for (m = 0; m < 5; ++m) {
                                if (first_5A_cluster[l] == second_5A_cluster[m]) break;
                            }
                            if (m == 5) {
                                if (k == 1) {
                                    k++;
                                    break;
                                }
                                unc[0] = first_5A_cluster[l];
                                k++;
                            }
                        }
                        if (k != 1) {
                            continue;
                        }

                        k = 0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                        for (l = 0; l < 3; ++l) {
                            if (second_5A_cluster[l] != cp[0] && second_5A_cluster[l] != cp[1]) {
                                for (m = 0; m < 5; ++m) {
                                    if (second_5A_cluster[l] == first_5A_cluster[m]) break;
                                }
                                if (m == 5) {
                                    if (k == 1) {
                                        k++;
                                        break;
                                    }
                                    unc[1] = second_5A_cluster[l];
                                    k++;
                                }
                            }
                        }
                        if (k != 1) {
                            continue;
                        }

                        for (int third_5A_pointer = second_5A_pointer + 1; third_5A_pointer < nmem_sp3c[first_5A_ring_particle]; ++third_5A_pointer) {
                            int third_5A_id = mem_sp3c[first_5A_ring_particle][third_5A_pointer];
                            if (third_5A_id > first_5A_id && third_5A_id > second_5A_id) {
                                int *third_5A_cluster = hcsp3c[third_5A_id];
                                n = 0;  // check common SP3 ring particles are exactly cp
                                for (l = 0; l < 3; ++l) {
                                    for (m = 0; m < 2; ++m) {
                                        if (cp[m] == third_5A_cluster[l]) {
                                            n++;
                                        }
                                    }
                                }
                                if (n != 2) continue;

                                n = 0;  // check spindles are exactly sother
                                for (l = 3; l < 5; ++l) {
                                    for (m = 0; m < 2; ++m) {
                                        if (sother[m] == third_5A_cluster[l]) {
                                            n++;
                                        }
                                    }
                                }
                                if (n != 2) continue;

                                n = 0;    // find particle from SP3 ring of mem_sp3c[sp3c[i][j2]][j] which is common with 5A sp3c_i
                                for (l = 0; l < 3; ++l) {
                                    if (third_5A_cluster[l] == cp[0] ||
                                        third_5A_cluster[l] == cp[1])
                                        continue;
                                    for (m = 0; m < 5; ++m) {
                                        if (third_5A_cluster[l] == first_5A_cluster[m] ||
                                            third_5A_cluster[l] ==
                                            second_5A_cluster[l])
                                            break;
                                    }
                                    if (m == 5) {
                                        if (n == 1) {
                                            n++;
                                            break;
                                        }
                                        unc[2] = third_5A_cluster[l];
                                        n++;
                                    }
                                }
                                if (n != 1) continue;

                                // Now we have found the 8K cluster
                                if (n8K == m8K) {
                                    hc8K = resize_2D_int(hc8K, m8K, m8K + incrStatic, clusSize, -1);
                                    m8K = m8K + incrStatic;
                                }
                                // hc8K key: (SP3_common_1, SP3_common_2, spindle_1, spindle_2, spindle_3, other_SP3_1, other_SP3_2, other_SP3_3)
                                hc8K[n8K][0] = cp[0];
                                hc8K[n8K][1] = cp[1];
                                hc8K[n8K][2] = scom;
                                hc8K[n8K][3] = sother[0];
                                hc8K[n8K][4] = sother[1];
                                hc8K[n8K][5] = unc[0];
                                hc8K[n8K][6] = unc[1];
                                hc8K[n8K][7] = unc[2];

                                quickSort(&hc8K[n8K][0], 2);
                                quickSort(&hc8K[n8K][2], 3);
                                quickSort(&hc8K[n8K][5], 3);

                                Cluster_Write_8K();
                            }
                        }
                    }
                }
            }
        }
    }
}

int is_particle_in_5A_ring(const int *second_5A_cluster, int first_5A_ring_particle, int *cp) {
    for (int ring_pointer = 0; ring_pointer < 3; ++ring_pointer) {
        if (first_5A_ring_particle == second_5A_cluster[ring_pointer]) {
            cp[0] = first_5A_ring_particle;
            return 1;
        }
    }
    return 0;
}

void Cluster_Write_8K() {
    if (s8K[hc8K[n8K][2]] == 'C') s8K[hc8K[n8K][2]] = 'B';
    if (s8K[hc8K[n8K][3]] == 'C') s8K[hc8K[n8K][3]] = 'B';
    if (s8K[hc8K[n8K][4]] == 'C') s8K[hc8K[n8K][4]] = 'B';
    if (s8K[hc8K[n8K][5]] == 'C') s8K[hc8K[n8K][5]] = 'B';
    if (s8K[hc8K[n8K][6]] == 'C') s8K[hc8K[n8K][6]] = 'B';
    if (s8K[hc8K[n8K][7]] == 'C') s8K[hc8K[n8K][7]] = 'B';
    s8K[hc8K[n8K][0]] = 'O';
    s8K[hc8K[n8K][1]] = 'O';
    ++n8K;
}