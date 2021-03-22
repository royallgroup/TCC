/* TCC: A topological cluster classification code with temporal tracking of clusters. */

#include "globals.h"
#include "bonds.h"
#include "rings.h"
#include "output.h"
#include "stats.h"
#include "tools.h"
#include "input.h"
#include "6Z.h"
#include "7K.h"
#include "7T.h"
#include "8A.h"
#include "8B.h"
#include "8K.h"
#include "9A.h"
#include "9B.h"
#include "9K.h"
#include "10A.h"
#include "10K.h"
#include "10W.h"
#include "11A.h"
#include "11B.h"
#include "11C.h"
#include "11W.h"
#include "11F.h"
#include "12A.h"
#include "12B.h"
#include "12E.h"
#include "12K.h"
#include "13B.h"
#include "HCP.h"
#include "FCC.h"
#include "BCC9.h"

char* cluster_names[] = {"sp3a", "sp3b", "sp3c", "sp4a", "sp4b", "sp4c", "sp5a", "sp5b", "sp5c",
                         "6A", "6Z", "7K", "7T_a", "7T_s", "8A", "8B", "8K", "9A", "9B", "9K", "10A", "10B",
                         "10K", "10W", "11A", "11B", "11C", "11E", "11F", "11W", "12A", "12B", "12D",
                         "12E", "12K", "13A", "13B", "13K", "FCC", "HCP", "BCC_9", "-1"};

int cluster_size[] = {3, 4, 5, 4, 5, 6, 5, 6, 7,
                      6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10,
                      10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12,
                      12, 12, 13, 13, 13, 13, 13, 9, -1};

int* do_cluster_list[] = {&dosp3a, &dosp3b, &dosp3c, &dosp4a, &dosp4b, &dosp4c, &dosp5a, &dosp5b, &dosp5c,
                          &do6A, &do6Z, &do7K, &do7T_a, &do7T_s, &do8A, &do8B, &do8K, &do9A, &do9B, &do9K, &do10A, &do10B,
                          &do10K, &do10W, &do11A, &do11B, &do11C, &do11E, &do11F, &do11W, &do12A, &do12B, &do12D,
                          &do12E, &do12K, &do13A, &do13B, &do13K, &doFCC, &doHCP, &doBCC9, NULL};

int* num_cluster_list[] = {&nsp3a, &nsp3b, &nsp3c, &nsp4a, &nsp4b, &nsp4c, &nsp5a, &nsp5b, &nsp5c,
                           &n6A, &n6Z, &n7K, &n7T_a, &n7T_s, &n8A, &n8B, &n8K, &n9A, &n9B, &n9K, &n10A, &n10B,
                           &n10K, &n10W, &n11A, &n11B, &n11C, &n11E, &n11F, &n11W, &n12A, &n12B, &n12D,
                           &n12E, &n12K, &n13A, &n13B, &n13K, &nFCC, &nHCP, &nBCC_9, NULL};

char** raw_list[] = {&ssp3a, &ssp3b, &ssp3c, &ssp4a, &ssp4b, &ssp4c, &ssp5a, &ssp5b, &ssp5c,
                     &s6A, &s6Z, &s7K, &s7T_a, &s7T_s, &s8A, &s8B, &s8K, &s9A, &s9B, &s9K, &s10A, &s10B,
                     &s10K, &s10W, &s11A, &s11B, &s11C, &s11E, &s11F, &s11W, &s12A, &s12B, &s12D,
                     &s12E, &s12K, &s13A, &s13B, &s13K, &sFCC, &sHCP, &sBCC_9, NULL};

int*** cluster_list[] = {&hcsp3a, &hcsp3b, &hcsp3c, &hcsp4a, &hcsp4b, &hcsp4c, &hcsp5a, &hcsp5b, &hcsp5c,
                         &hc6A, &hc6Z, &hc7K, &hc7T_a, &hc7T_s, &hc8A, &hc8B, &hc8K, &hc9A, &hc9B, &hc9K, &hc10A, &hc10B,
                         &hc10K, &hc10W, &hc11A, &hc11B, &hc11C, &hc11E, &hc11F, &hc11W, &hc12A, &hc12B, &hc12D,
                         &hc12E, &hc12K, &hc13A, &hc13B, &hc13K, &hcFCC, &hcHCP, &hcBCC_9, NULL};

int* cluster_list_width[] = {&msp3a, &msp3b, &msp3c, &msp4a, &msp4b, &msp4c, &msp5a, &msp5b, &msp5c,
                             &m6A, &m6Z, &m7K, &m7T_a, &m7T_s, &m8A, &m8B, &m8K, &m9A, &m9B, &m9K, &m10A, &m10B,
                             &m10K, &m10W, &m11A, &m11B, &m11C, &m11E, &m11F, &m11W, &m12A, &m12B, &m12D,
                             &m12E, &m12K, &m13A, &m13B, &m13K, &mFCC, &mHCP, &mBCC_9, NULL};

int main(int argc, char **argv) {
    int current_frame_number;
    struct xyz_info input_xyz_info;

    print_version_number();
    validate_cluster_lists();
    read_ini_file();
    read_clusters_to_analyse();
    analyse_cluster_dependencies();

    input_xyz_info = parse_xyz_file();
    frames_to_analyse = check_frame_numbers(input_xyz_info.total_frames);
    parse_box_file(frames_to_analyse);

    initialise_run_variables();
    Setup_Output_Files();

    for (current_frame_number = 0; current_frame_number < frames_to_analyse; current_frame_number++) {

        printf("Processing frame %d.\n", current_frame_number);
        particles_in_current_frame = input_xyz_info.num_particles[current_frame_number];
        initialise_frame_variables();
        get_box_size(current_frame_number);
        get_frame_coordinates_from_xyz(&input_xyz_info, current_frame_number);
        build_bond_network(current_frame_number);
        if (dosp3 == 1) get_basic_clusters();
        if (do6Z == 1) Clusters_Get6Z();
        if (do7K == 1) Clusters_Get7K();
        if (do7T_a == 1 || do7T_s == 1) Clusters_Get7T();
        if (do8A == 1) Clusters_Get8A();
        if (do8B == 1) Clusters_Get8B();
        if (do8K == 1) Clusters_Get8K();
        if (do9A == 1) Clusters_Get9A();
        if (do9B == 1) Clusters_Get9B_10B_11E_12D();
        if (do9K == 1) Clusters_Get9K();
        if (do10A == 1) Clusters_Get10A();
        if (do10K == 1) Clusters_Get10K();
        if (do10W == 1) Clusters_Get10W();
        if (do11A == 1) Clusters_Get11A();
        if (do11B == 1) Clusters_Get11B();
        if (do11C == 1) Clusters_Get11C();
        if (do11F == 1) Clusters_Get11F_13K();
        if (do11W == 1) Clusters_Get11W();
        if (do12A == 1) Clusters_Get12A();
        if (do12B == 1) Clusters_Get12B_13A();
        if (do12E == 1) Clusters_Get12E();
        if (do12K == 1) Clusters_Get12K();
        if (do13B == 1) Clusters_Get13B();
        if (doFCC == 1) Clusters_GetFCC();
        if (doHCP == 1) Clusters_GetHCP();
        if (doBCC9 == 1) Clusters_GetBCC_9();
        write_output_files(current_frame_number);
        free_frame_variables();
        printf("Cluster analysis for frame %d complete\n\n", current_frame_number);
    }

    // Do post analysis statistics
    count_mean_pop_per_frame(frames_to_analyse);
    Stats_Report();
    Write_Pop_Per_Frame(frames_to_analyse);

    free_xyz_info(&input_xyz_info);
    free_run_variables();

    printf("\n\nTCC completed successfully.\n\n");
    return 0;
}
