/* Alex Malins - alex.malins@gmail.com */
/* TCC: A topological cluster classification code with temporal tracking of clusters. */
/* Not for general consumption */

#include "globals.h"
#include "bonds.h"
#include "rings.h"
#include "output.h"
#include "stats.h"
#include "tools.h"
#include "input.h"

#include "6MW.h"
#include "6Z.h"

#include "7K.h"
#include "7MW.h"
#include "7PAB.h"
#include "7T.h"

#include "8A.h"
#include "8MW.h"
#include "8B.h"
#include "8K.h"
#include "8O.h"
#include "8PAA.h"
#include "8PAB.h"
#include "8PBB.h"

#include "9A.h"
#include "9B.h"
#include "9K.h"
#include "9MW.h"
#include "9PAA.h"
#include "9PAB.h"
#include "9PBB.h"
#include "9S.h"

#include "10A.h"
#include "10K.h"
#include "10MW.h"
#include "10W.h"

#include "11A.h"
#include "11B.h"
#include "11C.h"
#include "11E.h"
#include "11F.h"
#include "11MW.h"
#include "11W.h"
#include "11O.h"
#include "11PAA.h"
#include "11PAB.h"
#include "11PBB.h"
#include "11S.h"
#include "11SB.h"
#include "11W.h"

#include "12A.h"
#include "12B.h"
#include "12E.h"
#include "12K.h"
#include "12MW.h"
#include "12O.h"
#include "12PAA.h"
#include "12PAB.h"
#include "12PBB.h"
#include "12S.h"
#include "12SB.h"

#include "13B.h"
#include "13K.h"
#include "13MW.h"
#include "13PAB.h"
#include "13PBB.h"
#include "13S.h"
#include "13SB.h"

#include "14O.h"
#include "HCP.h"
#include "FCC.h"
#include "BCC9.h"

char* cluster_names[] = {"sp3a",  "sp3b",  "sp3c",  "sp4a", "sp4b", "sp4c",  "sp5a",  "sp5b", "sp5c",         // 9
                         "6A",    "6MW",   "6Z",                                                              // 12
                         "7K",    "7MW",   "7PAB",  "7T_a", "7T_s",                                           // 17
                         "8A",    "8B",    "8K",    "8MW",  "8O",   "8PAA",  "8PAB",  "8PBB",                 // 25
                         "9A",    "9B",    "9K",    "9MW",  "9PAA", "9PAB",  "9PBB",  "9S",                   // 33
                         "10A",   "10B",   "10K",   "10MW", "10O",  "10PAA", "10PAB", "10PBB", "10S", "10W",  // 43
                         "11A",   "11B",   "11C",   "11E",  "11F",  "11MW",  "11O",   "11PAA", "11PAB",       // 52
                         "11PBB", "11S",   "11SB",  "11W",                                                    // 56
                         "12A",   "12B",   "12D",   "12E",  "12K",  "12MW", "12O",                            // 63
                         "12PAA", "12PAB", "12PBB", "12S",  "12SB",                                           // 68
                         "13A",   "13B",   "13K",   "13MW",                                                   // 72
                         "13PAA", "13PAB", "13PBB", "13S",  "13SB",                                           // 77
                         "14O",   "FCC",   "HCP",   "BCC_9", "-1"};                                           // 81

int cluster_size[] = {3,  4,  5,  4,  5,  6,  5,  6,  7,        // 9
                      6,  6,  6,                                // 12
                      7,  7,  7,  7,  7,                        // 17
                      8,  8,  8,  8,  8,  8,  8,  8,            // 25
                      9,  9,  9,  9,  9,  9,  9,  9,            // 33
                      10, 10, 10, 10, 10, 10, 10, 10, 10, 10,   // 43
                      11, 11, 11, 11, 11, 11, 11, 11, 11,       // 52
                      11, 11, 11, 11,                           // 56
                      12, 12, 12, 12, 12, 12, 12,               // 63
                      12, 12, 12, 12, 12,                       // 68
                      13, 13, 13, 13,                           // 72
                      13, 13, 13, 13, 13,                       // 77
                      14, 13, 13, 9, -1};                       // 81

int* do_cluster_list[] = {&dosp3a,  &dosp3b,  &dosp3c,   &dosp4a,  &dosp4b, &dosp4c, &dosp5a,  &dosp5b, &dosp5c,            // 9
                          &do6A,    &do6MW,   &do6Z,                                                                        // 12
                          &do7K,    &do7MW,   &do7PAB,   &do7T_a,  &do7T_s,                                                 // 17
                          &do8A,    &do8B,    &do8K,     &do8MW,   &do8O,   &do8PAA, &do8PAB,  &do8PBB,                     // 25
                          &do9A,    &do9B,    &do9K,     &do9MW,   &do9PAA, &do9PAB, &do9PBB,  &do9S,                       // 33
                          &do10A,   &do10B,   &do10K,    &do10MW,  &do10O, &do10PAA, &do10PAB, &do10PBB, &do10S,   &do10W,  // 43
                          &do11A,   &do11B,   &do11C,    &do11E,   &do11F,  &do11MW, &do11O,   &do11PAA, &do11PAB,          // 52
                          &do11PBB, &do11S,   &do11SB,   &do11W,                                                            // 56
                          &do12A,   &do12B,   &do12D,    &do12E,   &do12K,  &do12MW, &do12O,                                // 63
                          &do12PAA, &do12PAB, &do12PBB,  &do12S,   &do12SB,                                                 // 68
                          &do13A,   &do13B,   &do13K,    &do13MW,                                                           // 72
                          &do13PAA, &do13PAB, &do13PBB,  &do13S,   &do13SB,                                                 // 77
                          &do14O,   &doFCC,   &doHCP,   &doBCC9,   NULL};                                                   // 81


int* num_cluster_list[] = {&nsp3a,   &nsp3b,  &nsp3c,  &nsp4a,  &nsp4b, &nsp4c, &nsp5a,  &nsp5b, &nsp5c,                   // 9
                           &n6A,     &n6MW,   &n6Z,                                                                        // 12
                           &n7K,     &n7MW,   &n7PAB,  &n7T_a,  &n7T_s,                                                    // 17
                           &n8A,     &n8B,    &n8K,    &n8MW,   &n8O,   &n8PAA, &n8PAB,  &n8PBB,                           // 25
                           &n9A,     &n9B,    &n9K,    &n9MW,   &n9PAA, &n9PAB, &n9PBB,  &n9S,                             // 33
                           &n10A,    &n10B,   &n10K,   &n10MW,  &n10O, &n10PAA, &n10PAB, &n10PBB, &n10S,   &n10W,          // 43
                           &n11A,    &n11B,   &n11C,   &n11E,   &n11F,  &n11MW, &n11O,   &n11PAA, &n11PAB,                 // 52
                           &n11PBB,  &n11S,   &n11SB,  &n11W,                                                              // 56
                           &n12A,    &n12B,   &n12D,   &n12E,   &n12K,  &n12MW, &n12O,                                     // 63
                           &n12PAA,  &n12PAB, &n12PBB, &n12S,   &n12SB,                                                    // 68
                           &n13A,    &n13B,   &n13K,   &n13MW,                                                             // 72
                           &n13PAA,  &n13PAB, &n13PBB, &n13S,   &n13SB,                                                    // 77
                           &n14O,    &nFCC,   &nHCP,   &nBCC_9, NULL};                                                     // 81

char** raw_list[] = {&ssp3a,  &ssp3b,  &ssp3c,  &ssp4a,   &ssp4b, &ssp4c,  &ssp5a, &ssp5b, &ssp5c,                // 9
                     &s6A,    &s6MW,   &s6Z,                                                                      // 12
                     &s7K,    &s7MW,   &s7PAB,  &s7T_a,   &s7T_s,                                                 // 17
                     &s8A,    &s8B,    &s8K,    &s8MW,    &s8O,   &s8PAA,  &s8PAB,  &s8PBB,                       // 25
                     &s9A,    &s9B,    &s9K,    &s9MW,    &s9PAA, &s9PAB,  &s9PBB,  &s9S,                         // 33
                     &s10A,   &s10B,   &s10K,   &s10MW,   &s10O,  &s10PAA, &s10PAB, &s10PBB, &s10S,    &s10W,     // 43
                     &s11A,   &s11B,   &s11C,   &s11E,    &s11F,  &s11MW,  &s11O,   &s11PAA, &s11PAB,             // 52
                     &s11PBB, &s11S,   &s11SB,  &s11W,                                                            // 56
                     &s12A,   &s12B,   &s12D,   &s12E,    &s12K,  &s12MW, &s12O,                                  // 63
                     &s12PAA, &s12PAB, &s12PBB, &s12S,    &s12SB,                                                 // 68
                     &s13A,   &s13B,   &s13K,   &s13MW,                                                           // 72
                     &s13PAA, &s13PAB, &s13PBB, &s13S,    &s13SB,                                                 // 77
                     &s14O,    &sFCC,   &sHCP,   &sBCC_9, NULL};                                                  // 81

int*** cluster_list[] = {&hcsp3a,  &hcsp3b,  &hcsp3c,  &hcsp4a, &hcsp4b, &hcsp4c, &hcsp5a, &hcsp5b, &hcsp5c,              // 9
                         &hc6A,    &hc6MW,   &hc6Z,                                                                       // 12
                         &hc7K,    &hc7MW,   &hc7PAB,  &hc7T_a, &hc7T_s,                                                  // 17
                         &hc8A,    &hc8B,    &hc8K,    &hc8MW,  &hc8O,   &hc8PAA,  &hc8PAB,  &hc8PBB,                     // 25
                         &hc9A,    &hc9B,    &hc9K,    &hc9MW,  &hc9PAA, &hc9PAB,  &hc9PBB,  &hc9S,                       // 33
                         &hc10A,   &hc10B,   &hc10K,   &hc10MW, &hc10O,  &hc10PAA, &hc10PAB, &hc10PBB, &hc10S,   &hc10W,  // 43
                         &hc11A,   &hc11B,   &hc11C,   &hc11E,  &hc11F,  &hc11MW,  &hc11O,   &hc11PAA, &hc11PAB,          // 52
                         &hc11PBB, &hc11S,   &hc11SB,  &hc11W,                                                            // 56
                         &hc12A,   &hc12B,   &hc12D,   &hc12E,  &hc12K,  &hc12MW,  &hc12O,                                // 63
                         &hc12PAA, &hc12PAB, &hc12PBB, &hc12S, &hc12SB,                                                   // 68
                         &hc13A,   &hc13B,   &hc13K,   &hc13MW,                                                           // 72
                         &hc13PAA, &hc13PAB, &hc13PBB, &hc13S, &hc13SB,                                                   // 77
                         &hc14O,   &hcFCC,   &hcHCP,   &hcBCC_9, NULL};                                                   // 81

int* cluster_list_width[] = {&msp3a,  &msp3b,  &msp3c,  &msp4a, &msp4b, &msp4c,  &msp5a,  &msp5b, &msp5c,           // 9
                             &m6A,    &m6MW,   &m6Z,                                                                // 12
                             &m7K,    &m7MW,   &m7PAB,  &m7T_a, &m7T_s,                                             // 17
                             &m8A,    &m8B,    &m8K,    &m8MW,  &m8O,   &m8PAA,  &m8PAB,  &m8PBB,                   // 25
                             &m9A,    &m9B,    &m9K,    &m9MW,  &m9PAA, &m9PAB,  &m9PBB,  &m9S,                     // 33
                             &m10A,   &m10B,   &m10K,   &m10MW, &m10O,  &m10PAA, &m10PAB, &m10PBB, &m10S,   &m10W,  // 43
                             &m11A,   &m11B,   &m11C,   &m11E,  &m11F,  &m11MW,  &m11O,   &m11PAA, &m11PAB,         // 52
                             &m11PBB, &m11S,   &m11SB,  &m11W,                                                      // 56
                             &m12A,   &m12B,   &m12D,   &m12E,  &m12K,  &m12MW,  &m12O,                             // 63
                             &m12PAA, &m12PAB, &m12PBB, &m12S,  &m12SB,                                             // 68
                             &m13A,   &m13B,   &m13K,   &m13MW,                                                     // 72
                             &m13PAA, &m13PAB, &m13PBB, &m13S,  &m13SB,                                             // 77
                             &m14O,   &mFCC,   &mHCP,   &mBCC_9, NULL};                                             // 81

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
        if (dosp3   == 1) get_basic_clusters();

        if (do6MW   == 1) Clusters_Get6MW();
        if (do6Z    == 1) Clusters_Get6Z();

        if (do7K    == 1) Clusters_Get7K();
        if (do7MW   == 1) Clusters_Get7MW();
        if (do7PAB  == 1) Clusters_Get7PAB();
        if (do7T_a  == 1 || do7T_s == 1) Clusters_Get7T();

        if (do8A    == 1) Clusters_Get8A();
        if (do8B    == 1) Clusters_Get8B();
        if (do8K    == 1) Clusters_Get8K();
        if (do8MW   == 1) Clusters_Get8MW();
        if (do8O    == 1) Clusters_Get8O();
        if (do8PAA  == 1) Clusters_Get8PAA();
        if (do8PBB  == 1) Clusters_Get8PBB();
        if (do8PAB  == 1) Clusters_Get8PAB();

        if (do9A    == 1) Clusters_Get9A();
        if (do9B    == 1) Clusters_Get9B_10B_11E_12D();
        if (do9K    == 1) Clusters_Get9K();
        if (do9MW   == 1) Clusters_Get9MW();

        if (do10A   == 1) Clusters_Get10A();
        if (do10K   == 1) Clusters_Get10K();
        if (do10MW  == 1) Clusters_Get10MW();
        if (do10W   == 1) Clusters_Get10W();

        if (do11A   == 1) Clusters_Get11A();
        if (do11B   == 1) Clusters_Get11B();
        if (do11C   == 1) Clusters_Get11C();
        if (do11F   == 1) Clusters_Get11F_13K();
        if (do11W   == 1) Clusters_Get11W();
        if (do11MW  == 1) Clusters_Get11MW();
        if (do11O   == 1) Clusters_Get11O();
        //        if (do11PAA == 1) Clusters_Get12PAA();
        //        if (do11PAB == 1) Clusters_Get12PAB();
        //        if (do11PBB == 1) Clusters_Get12PBB();
        if (do11S   == 1) Clusters_Get12S();
        if (do11SB  == 1) Clusters_Get12SB();

        if (do12A   == 1) Clusters_Get12A();
        if (do12B   == 1) Clusters_Get12B_13A();
        if (do12E   == 1) Clusters_Get12E();
        if (do12K   == 1) Clusters_Get12K();
        if (do12MW  == 1) Clusters_Get12MW();
        if (do12O   == 1) Clusters_Get12O();
        if (do12PAA == 1) Clusters_Get12PAA();
        if (do12PAB == 1) Clusters_Get12PAB();
        if (do12PBB == 1) Clusters_Get12PBB();
        if (do12S   == 1) Clusters_Get12S();
        if (do12SB  == 1) Clusters_Get12SB();

        if (do13A   == 1) Clusters_Get12B_13A();
        if (do13B   == 1) Clusters_Get13B();
        if (do13K   == 1) Clusters_Get11F_13K();
        if (do13MW  == 1) Clusters_Get13MW();
        if (do13PAA == 1) Clusters_Get13PAA();
        if (do13PAB == 1) Clusters_Get13PAB();
        if (do13PBB == 1) Clusters_Get13PBB();
        if (do13S   == 1) Clusters_Get13S();
        if (do13SB  == 1) Clusters_Get13SB();


        if (doFCC   == 1) Clusters_GetFCC();
        if (doHCP   == 1) Clusters_GetHCP();
        if (doBCC9  == 1) Clusters_GetBCC_9();
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
