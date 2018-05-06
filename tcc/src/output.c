#include <math.h>
#include <stats.h>
#include "output.h"
#include "tools.h"
#include "globals.h"

////////// Raw Writing //////////

void Write_Raw(int f) {
    int cluster_type;
    char file_name[200];
    FILE *file_pointer;

    for(cluster_type=0; cluster_type<num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            sprintf(file_name, "raw_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_%s",
                    fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[cluster_type]);
            file_pointer = fopen(file_name, "a");
            Write_Raw_Particle_Types(f, file_pointer, raw_list[cluster_type][0]);
            fclose(file_pointer);
        }
    }
}

void Write_Raw_Particle_Types(int f, FILE *thefile, const char *sarr) {
    int i;

    fprintf(thefile,"%ld\nframe %d\n", particles_in_current_frame, f + 1);
    for(i = 0; i < particles_in_current_frame; i++) {
        if (sarr[i] != 'C') {
            if (particle_type[i] == 1) {
                fprintf(thefile,"C\n");
            }
            else {
                fprintf(thefile,"D\n");
            }
        }
        else if (sarr[i] == 'C') {
            if (particle_type[i]==1){
                fprintf(thefile,"A\n");
            }
            else {
                fprintf(thefile,"B\n");
            }
        }
    }
}

////////// Bonds Writing //////////

void Write_Bonds_File(int f) {
    int i, j, sum;
    char errMsg[100];
    char output_file[200];
    FILE *bondsout;

    sum=0;
    for (i=0; i<particles_in_current_frame; ++i) {
        sum+=num_bonds[i];
    }
    if (sum%2!=0) {
        sprintf(errMsg,"Write_Bonds_File(): total number of bonds is not even %d\n",sum);
        exit(1);
    }

    sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,use_voronoi_bonds,fc,PBCs);
    bondsout=fopen(output_file, "a");

    fprintf(bondsout,"frame %d  total bonds %d\n",f,sum/2);
    for (i=0; i<particles_in_current_frame; ++i) {
        fprintf(bondsout,"%d    %d",i,num_bonds[i]);
        for (j=0; j<num_bonds[i]; ++j) {
            fprintf(bondsout,"  %d  %.5lg",bond_list[i][j],sqrt(squared_bondlengths[i][j]));
        }
        fprintf(bondsout,"\n");
    }
    fclose(bondsout);
}

////////// Centers Writing //////////

void Write_Cluster_Centers_xyz(int f, int cluster_type) {

    int num_centers = 0;
    FILE *output_file;
    char file_name[200];

    sprintf(file_name, "centers_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.%s_cen.xyz",
            fXmolName,rcutAA,rcutAB,rcutBB,use_voronoi_bonds,fc,PBCs, cluster_names[cluster_type]);

    output_file = fopen(file_name, "a");

    for(int particle_number = 0; particle_number < particles_in_current_frame; particle_number++) {
        if((*raw_list[cluster_type])[particle_number]=='S') {
            num_centers++;
        }
    }

    fprintf(output_file, "%d\nframe %d of %d\n", num_centers, f+1, frames_to_analyse);
    for(int particle_number = 0; particle_number < particles_in_current_frame; particle_number++) {
        if ((*raw_list[cluster_type])[particle_number] == 'S') {
            fprintf(output_file ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[particle_number], y[particle_number], z[particle_number]);
        }
    }

    fclose(output_file);
}

////////// Cluster XYZ Writing //////////

void Write_Cluster_XYZ(int f) {
    int i, cluster_type, num_particles;
    char output_file[200];
    FILE *file_pointer;

    for(cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            sprintf(output_file, "cluster_xyzs/%s.%s_clusts.xyz", fXmolName, cluster_names[cluster_type]);
            file_pointer = open_file(output_file, "a");

            num_particles = 0;
            for(i = 0; i < particles_in_current_frame; i++) {
                if((*raw_list[cluster_type])[i] != 'C') {
                    num_particles++;
                }
            }

            fprintf(file_pointer,"%d\n", num_particles);
            fprintf(file_pointer,"Frame number %d\n", f);
            for(i = 0; i < particles_in_current_frame; i++) {
                if ((*raw_list[cluster_type])[i] != 'C') {
                    fprintf(file_pointer, "%c\t%f\t%f\t%f\n", (*raw_list[cluster_type])[i], x[i], y[i], z[i]);
                }
            }
        }
    }
}

////////// Cluster ids Writing //////////

void Write_Cluster(int f) {
    int cluster_type;

    for(cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            Write_Cluster_Compostions(f, cluster_type);
        }
    }
}

void Write_Cluster_Compostions(int f, int cluster_type) {
    int i, j;
    char output_file[200];
    FILE *file_pointer;
    int num_clusters, clusSize;
    int **hc;

    num_clusters = *num_cluster_list[cluster_type];
    hc = *cluster_list[cluster_type];
    clusSize = cluster_size[cluster_type];

    sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
            fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[cluster_type]);
    file_pointer = open_file(output_file, "a");
    num_sort_columns = *num_cluster_list[cluster_type];
    qsort(*cluster_list[cluster_type], (size_t)num_sort_columns, sizeof(int *), sort_list_of_lists_of_ints);
    fprintf(file_pointer,"Frame Number %d\n",f);
    for (i = 0; i < num_clusters; i++) {
        for (j = 0; j < clusSize - 1; j++) {
            fprintf(file_pointer,"%d\t",hc[i][j]);
        }
        fprintf(file_pointer,"%d\n",hc[i][clusSize - 1]);
    }
    fclose(file_pointer);
}

////////// Pop per frame writing //////////

void Write_Pop_Per_Frame(int f) {
    char errMsg[1000], output[1000];
    int cluster_type, frame;
    FILE *file_pointer;

    if (doWritePopPerFrame==1) {

        sprintf(output, "%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame", fXmolName, rcutAA, rcutAB,
                rcutBB, use_voronoi_bonds, fc, PBCs);
        file_pointer = open_file(output, "a");

        if (file_pointer == NULL) {
            sprintf(errMsg, "main() : Error opening file %s", output);    // Always test file open
            Error(errMsg);
        }

        // Print column headings
        fprintf(file_pointer, "frame\t");
        for (frame = 0; frame < f; frame++) {
            fprintf(file_pointer, "%d\t", frame);
        }
        fprintf(file_pointer, "\n");

        // Print populations of each cluster type, one cluster per line
        for (cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {
            fprintf(file_pointer, "%s\t", cluster_names[cluster_type]);
            for (frame = 0; frame < f; frame++) {
                if (*do_cluster_list[cluster_type] == 1) {
                    fprintf(file_pointer, "%.15lg	", pop_per_frame[cluster_type][frame]);
                }
                else {
                    fprintf(file_pointer, "NA	");
                }
            }
            fprintf(file_pointer, "\n");
        }
        fclose(file_pointer);
    }
}

void write_output_files(int current_frame_number) {
    count_number_of_clusters();
    count_frame_cluster_population(current_frame_number);

    if (doWriteClus == 1) Write_Cluster(current_frame_number);
    if (doWriteRaw == 1) Write_Raw(current_frame_number);
    if (doWriteXYZ == 1) Write_Cluster_XYZ(current_frame_number);
    if (do11AcenXyz == 1) Write_Cluster_Centers_xyz(current_frame_number, eleven_A_number);
    if (do13AcenXyz == 1) Write_Cluster_Centers_xyz(current_frame_number, thirteen_A_number);
}