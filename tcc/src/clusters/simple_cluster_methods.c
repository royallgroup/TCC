#include "simple_cluster_methods.h"

int count_common_spindles_between_5As(const int *first_5A_cluster, const int *second_5A_cluster, int *scom) {
    int num_common_spindles = 0;
    for (int ring_1_pointer = 3; ring_1_pointer < 5; ring_1_pointer++) {
        for (int ring_2_pointer = 3; ring_2_pointer < 5; ring_2_pointer++) {
            if (first_5A_cluster[ring_1_pointer] == second_5A_cluster[ring_2_pointer]) {
                *scom = first_5A_cluster[ring_1_pointer];
                num_common_spindles++;
            }
        }
    }
    return num_common_spindles;
}