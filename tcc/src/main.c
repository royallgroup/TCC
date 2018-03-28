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
#include "6Z.h"
#include "7K.h"
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
#include "BCC15.h"
#include "BCC9.h"

int num_cluster_types = 39;

char* cluster_names[] = {"sp3a", "sp3b", "sp3c", "sp4a", "sp4b", "sp4c", "sp5a", "sp5b", "sp5c",
                         "6Z", "7K", "8A", "8B", "8K", "9A", "9B", "9K", "10A", "10B", "10K", "10W",
                         "11A", "11B", "11C", "11E", "11F", "11W", "12A", "12B", "12D",
                         "12E", "12K", "13A", "13B", "13K", "FCC", "HCP", "BCC_9", "BCC_15"};

int cluster_size[] = {3, 4, 5, 4, 5, 6, 5, 6, 7,
                      6, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10,
                      11, 11, 11, 11, 11, 11, 12, 12, 12,
                      12, 12, 13, 13, 13, 13, 13, 9, 15};

int* do_cluster_list[] = {&dosp3a, &dosp3b, &dosp3c, &dosp4a, &dosp4b, &dosp4c, &dosp5a, &dosp5b, &dosp5c,
                          &do6Z, &do7K, &do8A, &do8B, &do8K, &do9A, &do9B, &do9K, &do10A, &do10B, &do10K, &do10W,
                          &do11A, &do11B, &do11C, &do11E, &do11F, &do11W, &do12A, &do12B, &do12D,
                          &do12E, &do12K, &do13A, &do13B, &do13K, &doFCC, &doHCP, &doBCC9, &doBCC15};

int* num_cluster_list[] = {&nsp3a, &nsp3b, &nsp3c, &nsp4a, &nsp4b, &nsp4c, &nsp5a, &nsp5b, &nsp5c,
                           &n6Z, &n7K, &n8A, &n8B, &n8K, &n9A, &n9B, &n9K, &n10A, &n10B, &n10K, &n10W,
                           &n11A, &n11B, &n11C, &n11E, &n11F, &n11W, &n12A, &n12B, &n12D,
                           &n12E, &n12K, &n13A, &n13B, &n13K, &nFCC, &nHCP, &nBCC_9, &nBCC_15};

char** raw_list[] = {&ssp3a, &ssp3b, &ssp3c, &ssp4a, &ssp4b, &ssp4c, &ssp5a, &ssp5b, &ssp5c,
                     &s6Z, &s7K, &s8A, &s8B, &s8K, &s9A, &s9B, &s9K, &s10A, &s10B, &s10K, &s10W,
                     &s11A, &s11B, &s11C, &s11E, &s11F, &s11W, &s12A, &s12B, &s12D,
                     &s12E, &s12K, &s13A, &s13B, &s13K, &sFCC, &sHCP, &sBCC_9, &sBCC_15};

int*** cluster_list[] = {&hcsp3a, &hcsp3b, &hcsp3c, &hcsp4a, &hcsp4b, &hcsp4c, &hcsp5a, &hcsp5b, &hcsp5c,
                         &hc6Z, &hc7K, &hc8A, &hc8B, &hc8K, &hc9A, &hc9B, &hc9K, &hc10A, &hc10B, &hc10K, &hc10W,
                         &hc11A, &hc11B, &hc11C, &hc11E, &hc11F, &hc11W, &hc12A, &hc12B, &hc12D,
                         &hc12E, &hc12K, &hc13A, &hc13B, &hc13K, &hcFCC, &hcHCP, &hcBCC_9, &hcBCC_15};

int main(int argc, char **argv) {
    int current_frame_number, f, i;
    int remainder;
    struct xyz_info input_xyz_info;
    char fInputParamsName[50];

    sprintf(fInputParamsName,"inputparameters.ini");
    Setup_ReadIniFile(fInputParamsName);    // read input params
    printf("box size file: %s\n",fBoxSizeName);

    input_xyz_info = parse_xyz_file();

    parse_box_file(input_xyz_info.total_frames);

    max_particle_number = get_max_particle_number(input_xyz_info);

    Initialise_Global_Variables();
    Stats_Init();
    Setup_Output_Files();

    f=0;
    for (current_frame_number=0; current_frame_number < input_xyz_info.total_frames; current_frame_number++) {
        remainder=current_frame_number%SAMPLEFREQ;

        if (remainder==0 && f<FRAMES) {
            current_frame_particle_number = input_xyz_info.num_particles[current_frame_number];
            Reset_Frame_Variables();
            if (box_type != 1) get_box_size(current_frame_number);
            get_xyz_frame(&input_xyz_info, current_frame_number);
            Get_Bonds();
            if (doWriteBonds == 1) Write_Bonds_File(f);

            for (i = 0; i < current_frame_particle_number; i++) {
                if (num_bonds[i] > maxnb) maxnb = num_bonds[i];
            }
            if (dosp3 == 1) get_basic_clusters();
            if (do6Z == 1) Clusters_Get6Z();
            if (do7K == 1) Clusters_Get7K();
            if (do8A == 1) Clusters_Get8A();
            if (do8B == 1) Clusters_Get8B();
            if (do8K == 1) Clusters_Get8K();
            if (do9A == 1) Clusters_Get9A();
            if (do9B == 1) Clusters_Get9B_10B_11B_11E_12D();
            if (do9K == 1) Clusters_Get9K();
            if (do10A == 1) Clusters_Get10A();
            if (do10K == 1) Clusters_Get10K();
            if (do10W == 1) Clusters_Get10W();
            if (do11A == 1) Clusters_Get11A();
            if (do11C == 1) Clusters_Get11C();
            if (do11F == 1) Clusters_Get11F_12E_13K();
            if (do11W == 1) Clusters_Get11W();
            if (do12A == 1) Clusters_Get12A();
            if (do12B == 1) Clusters_Get12B_13A();
            if (do12E == 1) Clusters_Get12E();
            if (do12K == 1) Clusters_Get12K();
            if (do13B == 1) Clusters_Get13B();
            if (doFCC == 1) Clusters_GetFCC();
            if (doHCP == 1) Clusters_GetHCP();
            if (doBCC9 == 1) Clusters_GetBCC_9();
            if (doBCC15 == 1) Clusters_GetBCC_15();

            // Write output files
            count_number_of_clusters();
            Pop_Per_Frame(f);

            if (doWriteClus == 1) Write_Cluster(f);
            if (doWriteRaw == 1) Write_Raw(f);
            if (doWriteXYZ == 1) Write_Cluster_XYZ(f);
            if (do11AcenXyz == 1) Write_Cluster_Centers_xyz(f, 21);
            if (do13AcenXyz == 1) Write_Cluster_Centers_xyz(f, 32);

            if (USELIST==1) {
                free(head);
                free(linked_list);
            }
            
            printf("f%d complete\n", f);
            f++;
            if (f == FRAMES) break;
        }
    }

    // Do post analysis statistics
    count_mean_pop_per_frame(f);
    Stats_Report();
    if (doWritePopPerFrame==1) {
        Write_Pop_Per_Frame(f);
    }

    free(input_xyz_info.num_particles);
    free(input_xyz_info.frame_offsets);
    Stats_FreeMem();
    Free_All_Variables();

    printf("\n\nFIN \n\n");
    return 0;
}