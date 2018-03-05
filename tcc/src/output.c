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
                    fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_type]);
            file_pointer = fopen(file_name, "a");
            Write_Raw_Particle_Types(f, file_pointer, raw_list[cluster_type][0]);
            fclose(file_pointer);
        }
    }
}

void Write_Raw_Particle_Types(int f, FILE *thefile, const char *sarr) {
    int i;

    fprintf(thefile,"%d\nframe %d\n",current_frame_particle_number,f+1);
    for(i=0; i<current_frame_particle_number; i++) {
        if (sarr[i]!='C') {
            if (particle_type[i]==1) fprintf(thefile,"C\n");
            else fprintf(thefile,"D\n");
        }
        else if (sarr[i]=='C') {
            if (particle_type[i]==1) fprintf(thefile,"A\n");
            else fprintf(thefile,"B\n");
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
    for (i=0; i<current_frame_particle_number; ++i) {
        sum+=num_bonds[i];
    }
    if (sum%2!=0) {
        sprintf(errMsg,"Write_Bonds_File(): total number of bonds is not even %d\n",sum);
        exit(1);
    }

    sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    bondsout=fopen(output_file, "a");

    fprintf(bondsout,"frame %d  total bonds %d\n",f,sum/2);
    for (i=0; i<current_frame_particle_number; ++i) {
        fprintf(bondsout,"%d    %d",i,num_bonds[i]);
        for (j=0; j<num_bonds[i]; ++j) {
            fprintf(bondsout,"  %d  %.5lg",bNums[i][j],bondlengths[i][j]);
        }
        fprintf(bondsout,"\n");
    }
    fclose(bondsout);
}

////////// Centers Writing //////////

void Write_Cluster_Centers_xyz(int f, int cluster_type) {

    int i, num_centers=0;
    FILE *output_file;
    char file_name[200];

    sprintf(file_name, "centers_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.%s_cen.xyz",
            fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, cluster_names[cluster_type]);

    output_file = fopen(file_name, "a");

    for(i=0; i<current_frame_particle_number; i++) if((*raw_list[cluster_type])[i]=='S') ++num_centers;

    fprintf(output_file,"%d\nframe %d of %d\n",num_centers,f+1,FRAMES);
    for(i=0; i<current_frame_particle_number; i++) {
        if ((*raw_list[cluster_type])[i]=='S') {
            fprintf(output_file ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
        }
    }

    fclose(output_file);
}

////////// Cluster Writing //////////

void Write_Cluster(int f) {
    int cluster_type;

    for(cluster_type=0; cluster_type<num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            Write_Cluster_Compostions(f, *num_cluster_list[cluster_type], *cluster_list[cluster_type],
                                      cluster_size[cluster_type], cluster_type);
        }
    }
}

void Write_Cluster_Compostions(int f, int num_clusters, int **hc, int clusSize, int cluster_number) {
    int i,j;
    char output_file[200];
    FILE *file_pointer;

    sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
            fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
    file_pointer = open_file(output_file, "a");

    fprintf(file_pointer,"Frame Number %d\n",f);
    for (i=0;i<num_clusters;i++) {
        fprintf(file_pointer,"%d",hc[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hc[i][j]);
        fprintf(file_pointer,"	%d\n",hc[i][clusSize-1]);
    }
    fclose(file_pointer);
}

////////// Pop per frame writing //////////

void Write_Pop_Per_Frame(int f) {
    char errMsg[1000], output[1000];
    int i;
    FILE *file_pointer;

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_pointer = open_file(output, "a");

    if (file_pointer==NULL)  {
        sprintf(errMsg,"main() : Error opening file %s",output);	// Always test file open
        Error(errMsg);
    }
    fprintf(file_pointer,"%s\n",output);

    fprintf(file_pointer,"frame	");
    for(i=0; i<num_cluster_types; i++) {
        fprintf(file_pointer, "%s	", cluster_names[i]);
    }
    fprintf(file_pointer,"\n");

    fprintf(file_pointer,"mean	");
    for(i=0; i<num_cluster_types; i++) {
        if(*do_cluster_list[i] == 1) {
            fprintf(file_pointer, "%.15lg	", mean_pop_per_frame[i]);
        }
        else {
            fprintf(file_pointer, "NA	");
        }
    }
    fprintf(file_pointer,"\n");

    for (f=0;f<FRAMES;f++) {
        fprintf(file_pointer,"%d	",f);
        for(i=0; i<num_cluster_types; i++) {
            if(*do_cluster_list[i] == 1) {
                fprintf(file_pointer, "%.15lg	", pop_per_frame[i][f]);
            }
            else {
                fprintf(file_pointer, "NA	");
            }
        }
        fprintf(file_pointer,"\n");
    }
    fclose(file_pointer);
}
