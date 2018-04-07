#include <globals.h>
#include <tools.h>
#include "7T.h"

int get_new_particle(const int *new_5A_cluster, int spindle_id);

int check_spindle_is_new_5A_spindle(int spindle_id, const int *new_5A_cluster);

void Clusters_Get7T() {
    int old_6Z_id, new_5A_id;
    int *old_6Z_cluster, *new_5A_cluster;
    int spindle_id, new_particle_id;
    int new_5A_pointer;
    int bond_counter;

    for (old_6Z_id = 0; old_6Z_id < n6Z; ++old_6Z_id) {
        old_6Z_cluster = hc6Z[old_6Z_id];
        for (int spindle_pointer = 4; spindle_pointer < 6; ++spindle_pointer) {
            spindle_id = old_6Z_cluster[spindle_pointer];
            for (new_5A_pointer = 0; new_5A_pointer < nmem_sp3c[spindle_id]; ++new_5A_pointer) {
                new_5A_id = mem_sp3c[spindle_id][new_5A_pointer];
                new_5A_cluster = hcsp3c[new_5A_id];

                if (check_spindle_is_new_5A_spindle(spindle_id, new_5A_cluster) == 1 ) {
                    bond_counter = check_ring_bonds(new_5A_cluster, old_6Z_cluster);
                    if (bond_counter == 21 || bond_counter == 22 || bond_counter == 25) {
                        new_particle_id = get_new_particle(new_5A_cluster, spindle_id);
                        check_7T_type(bond_counter, old_6Z_cluster, new_particle_id);
                    }
                }
            }
        }
    }
}

int check_spindle_is_new_5A_spindle(int spindle_id, const int *new_5A_cluster) {
    if (new_5A_cluster[3] == spindle_id) {
        return 1;
    }
    else if (new_5A_cluster[4] == spindle_id) {
        return 1;
    }
    else {
        return 0;
    }
}

int get_new_particle(const int *new_5A_cluster, int spindle_id) {
    if (new_5A_cluster[3] == spindle_id) return new_5A_cluster[4];
    else if (new_5A_cluster[4] == spindle_id) return new_5A_cluster[3];
    else Error("New spindle not found.");
}

int check_ring_bonds(const int *new_5A_cluster, const int *old_6Z_cluster) {
    int bond_counter = 0;

    for (int ring_pointer = 0; ring_pointer < 3; ++ring_pointer) {
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[0]) {
            bond_counter += 1;
            continue;
        }
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[1]) {
            bond_counter += 2;
            continue;
        }
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[2]) {
            bond_counter += 4;
            continue;
        }
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[3]) {
            bond_counter += 8;
            continue;
        }
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[4]) {
            bond_counter += 16;
            continue;
        }
        if (new_5A_cluster[ring_pointer] == old_6Z_cluster[5]) {
            bond_counter += 16;
            continue;
        }
    }
    return bond_counter;
}

void check_7T_type(int bond_counter, const int *old_6Z, int new_particle) {
    if (bond_counter == 21) {
        //detected symmetric
        add_7T_s(old_6Z, new_particle);
    }
    if (bond_counter == 22) {
        // detected assymetric type 1
        add_7T_a(old_6Z, new_particle);
    }
    if (bond_counter == 25) {
        // detected assymetric type 2
        add_7T_a(old_6Z, new_particle);
    }
}

void add_7T_a(const int *old_6Z, int new_particle) {
    int clusSize = 7;

    if (n7T_a == m7T_a) {
        hc7T_a = resize_2D_int(hc7T_a, m7T_a, m7T_a + incrStatic, clusSize, -1);
        m7T_a = m7T_a + incrStatic;
    }

    hc7T_a[n7T_a][0] = old_6Z[0];
    hc7T_a[n7T_a][1] = old_6Z[1];
    hc7T_a[n7T_a][2] = old_6Z[2];
    hc7T_a[n7T_a][3] = old_6Z[3];
    hc7T_a[n7T_a][4] = old_6Z[4];
    hc7T_a[n7T_a][5] = old_6Z[5];
    hc7T_a[n7T_a][6] = new_particle;

    for (int i = 0; i < 7; ++i) {
        s7T_a[hc7T_a[n7T_a][i]] = 'B';
    }
    ++n7T_a;
}

void add_7T_s(const int *old_6Z, int new_particle) {
    int clusSize = 7;

    if (n7T_s == m7T_s) {
        hc7T_s = resize_2D_int(hc7T_s, m7T_s, m7T_s + incrStatic, clusSize, -1);
        m7T_s = m7T_s
                + incrStatic;
    }

    hc7T_s[n7T_s][0] = old_6Z[0];
    hc7T_s[n7T_s][1] = old_6Z[1];
    hc7T_s[n7T_s][2] = old_6Z[2];
    hc7T_s[n7T_s][3] = old_6Z[3];
    hc7T_s[n7T_s][4] = old_6Z[4];
    hc7T_s[n7T_s][5] = old_6Z[5];
    hc7T_s[n7T_s][6] = new_particle;

    for (int i = 0; i < 7; ++i) {
        s7T_s[hc7T_s[n7T_s][i]] = 'B';
    }
    ++n7T_s;
}
