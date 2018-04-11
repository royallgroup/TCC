#include "stats.h"
#include "globals.h"
#include "tools.h"

void count_number_of_clusters() {
    for(int cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        total_clusters[cluster_type] += *num_cluster_list[cluster_type];
    }
}

void Stats_Report() {
    // The report which is printed to the console and static_clust.
    char output_name[200];
    FILE *output_file;

    sprintf(output_name, "%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.static_clust", fXmolName, rcutAA, rcutAB,
            rcutBB, use_voronoi_bonds, fc, PBCs);

    output_file = open_stats_report_file(output_name);

    stats_report_title(output_name, output_file);
    stats_report_clusters(output_file);


    fclose(output_file);
}

FILE *open_stats_report_file(const char *output_name) {
    char errMsg[1000];
    FILE *output_file;

    output_file = fopen(output_name,"w");
    if (output_file==NULL)  {
        sprintf(errMsg,"Stats_Report(): Error opening file %s",output_name);	// Always test file open
        Error(errMsg);
    }
    return output_file;
}

void stats_report_clusters(FILE *output_file) {
    char buffer[1000];

    for(int cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {
        if(*do_cluster_list[cluster_type] == 1) {
            sprintf(buffer, "%s	%d	%d	%.5lg\n", cluster_names[cluster_type], total_clusters[cluster_type],
                    num_gross_particles[cluster_type], mean_pop_per_frame[cluster_type]);
        }
        else {
            sprintf(buffer, "%s	NA	NA	NA\n", cluster_names[cluster_type]);
        }
        printf("%s", buffer);
        fprintf(output_file, "%s", buffer);
    }
}

void stats_report_title(const char *output_name, FILE *output_file) {
    char buffer[1000];

    fprintf(output_file, "%s\n", output_name);
    sprintf(buffer, "Cluster type	Number of clusters	Gross particles	Mean Pop Per Frame\n");
    printf("%s", buffer);
    fprintf(output_file, "%s", buffer);
}

void count_frame_cluster_population(int f) {

    for(int cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {
        pop_per_frame[cluster_type][f] = 0;
        for(int particle_number = 0; particle_number < particles_in_current_frame; particle_number++) {
            if ((*raw_list[cluster_type])[particle_number] != 'C') {
                pop_per_frame[cluster_type][f] += 1;
                num_gross_particles[cluster_type]++;
            }
        }

        // Normalise cluster count
        pop_per_frame[cluster_type][f] /= particles_in_current_frame;
    }
}

void count_mean_pop_per_frame(int frames_analysed) {

    for(int cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        mean_pop_per_frame[cluster_type] = 0;
        for (int frame_number = 0; frame_number < frames_analysed; frame_number++) {
            mean_pop_per_frame[cluster_type] += pop_per_frame[cluster_type][frame_number];
        }
        mean_pop_per_frame[cluster_type] /= frames_analysed;
    }
}
