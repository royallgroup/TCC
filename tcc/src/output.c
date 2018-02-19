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
            Write_Raw_Xmol(f, file_pointer, raw_cluster_list[cluster_type][0]);
            fclose(file_pointer);
        }
    }
}

void Write_Raw_Xmol(int f, FILE *thefile, const char *sarr) {
    int i;

    fprintf(thefile,"%d\nframe %d of %d\n",N,f+1,TOTALFRAMES);
    for(i=0; i<N; i++) {
        if (sarr[i]!='C') {
            if (rtype[i]==1) fprintf(thefile,"C\n");
            else fprintf(thefile,"D\n");
        }
        else if (sarr[i]=='C') {
            if (rtype[i]==1) fprintf(thefile,"A\n");
            else fprintf(thefile,"B\n");
        }
    }
}

////////// Bonds Writing //////////

void Bonds_WriteBonds(int f) {
    int i, j, sum;
    char errMsg[100];
    char output_file[200];
    FILE *bondsout;

    sum=0;
    for (i=0; i<N; ++i) {
        sum+=cnb[i];
    }
    if (sum%2!=0) {
        sprintf(errMsg,"Bonds_WriteBonds(): total number of bonds is not even %d\n",sum);
        exit(1);
    }

    sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    bondsout=fopen(output_file, "a");

    fprintf(bondsout,"frame %d  total bonds %d\n",f,sum/2);
    for (i=0; i<N; ++i) {
        fprintf(bondsout,"%d    %d",i,cnb[i]);
        for (j=0; j<cnb[i]; ++j) {
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

    for(i=0; i<N; i++) if((*raw_cluster_list[cluster_type])[i]=='S') ++num_centers;

    fprintf(output_file,"%d\nframe %d of %d\n",num_centers,f,FRAMES);
    for(i=0; i<N; i++) {
        if ((*raw_cluster_list[cluster_type])[i]=='S') {
            fprintf(output_file ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
        }
    }

    fclose(output_file);
}

////////// Cluster Writing //////////

void Write_Cluster_Xmol(int f, int num_clusters, int **hc, int clusSize, int cluster_number) {
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

void Write_Cluster_sp3(int f, int cluster_number) {
    int i,j;
    int clusSize=3;
    char output_file[200];
    FILE* file_pointer;

    sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
            fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
    file_pointer = open_file(output_file, "a");

    fprintf(file_pointer,"Frame Number%d\n",f);
    for (i=0;i<nsp3a;i++) {
        fprintf(file_pointer,"%d",hcsp3a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp3a[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp3a[i][clusSize-1]);
    }
    for (i=0;i<nsp3b;i++) {
        fprintf(file_pointer,"%d",hcsp3b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp3b[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp3b[i][clusSize-1]);
    }
    for (i=0;i<nsp3c;i++) {
        fprintf(file_pointer,"%d",hcsp3c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp3c[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp3c[i][clusSize-1]);
    }
    fclose(file_pointer);
}

void Write_Cluster_sp4(int f, int cluster_number) {
    int i,j;
    int clusSize=4;
    char output_file[200];
    FILE* file_pointer;

    sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
            fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
    file_pointer = open_file(output_file, "a");


    fprintf(file_pointer,"Frame Number%d\n",f);
    for (i=0;i<nsp4a;i++) {
        fprintf(file_pointer,"%d",hcsp4a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp4a[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp4a[i][clusSize-1]);
    }
    for (i=0;i<nsp4b;i++) {
        fprintf(file_pointer,"%d",hcsp4b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp4b[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp4b[i][clusSize-1]);
    }
    for (i=0;i<nsp4c;i++) {
        fprintf(file_pointer,"%d",hcsp4c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(file_pointer,"	%d",hcsp4c[i][j]);
        fprintf(file_pointer,"	%d\n",hcsp4c[i][clusSize-1]);
    }
}

void Write_Cluster_sp5(int f, int cluster_number) {
    int i, j;
    int clusSize = 5;
    char output_file[200];
    FILE* file_pointer;

    sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
            fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
    file_pointer = open_file(output_file, "a");

    fprintf(file_pointer,"Frame Number%d\n",f);
    for (i = 0; i < nsp5a; i++) {
        fprintf(file_pointer, "%d", hcsp5a[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(file_pointer, "	%d", hcsp5a[i][j]);
        fprintf(file_pointer, "	%d\n", hcsp5a[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5b; i++) {
        fprintf(file_pointer, "%d", hcsp5b[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(file_pointer, "	%d", hcsp5b[i][j]);
        fprintf(file_pointer, "	%d\n", hcsp5b[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5c; i++) {
        fprintf(file_pointer, "%d", hcsp5c[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(file_pointer, "	%d", hcsp5c[i][j]);
        fprintf(file_pointer, "	%d\n", hcsp5c[i][clusSize - 1]);
    }
}

void Write_Cluster(int f) {
    if(*do_cluster_list[0] == 1) Write_Cluster_sp3(f, 0);
    if(*do_cluster_list[1] == 1) Write_Cluster_Xmol(f, *num_cluster_list[1], hcsp3a, 3, 1);
    if(*do_cluster_list[2] == 1) Write_Cluster_Xmol(f, *num_cluster_list[2], hcsp3b, 4, 2);
    if(*do_cluster_list[3] == 1) Write_Cluster_Xmol(f, *num_cluster_list[3], hcsp3c, 5, 3);
    if(*do_cluster_list[4] == 1) Write_Cluster_sp4(f, 4);
    if(*do_cluster_list[5] == 1) Write_Cluster_Xmol(f, *num_cluster_list[5], hcsp4a, 4, 5);
    if(*do_cluster_list[6] == 1) Write_Cluster_Xmol(f, *num_cluster_list[6], hcsp4b, 5, 6);
    if(*do_cluster_list[7] == 1) Write_Cluster_Xmol(f, *num_cluster_list[7], hcsp4c, 6, 7);
    if(*do_cluster_list[8] == 1) Write_Cluster_sp5(f, 8);
    if(*do_cluster_list[9] == 1) Write_Cluster_Xmol(f, *num_cluster_list[9], hcsp5a, 5, 9);
    if(*do_cluster_list[10] == 1) Write_Cluster_Xmol(f, *num_cluster_list[10], hcsp5b, 6, 10);
    if(*do_cluster_list[11] == 1) Write_Cluster_Xmol(f, *num_cluster_list[11], hcsp5c, 7, 11);
    if(*do_cluster_list[12] == 1) Write_Cluster_Xmol(f, *num_cluster_list[12], hc6Z, 6, 12);
    if(*do_cluster_list[13] == 1) Write_Cluster_Xmol(f, *num_cluster_list[13], hc7K, 7, 13);
    if(*do_cluster_list[14] == 1) Write_Cluster_Xmol(f, *num_cluster_list[14], hc8A, 8, 14);
    if(*do_cluster_list[15] == 1) Write_Cluster_Xmol(f, *num_cluster_list[15], hc8B, 8, 15);
    if(*do_cluster_list[16] == 1) Write_Cluster_Xmol(f, *num_cluster_list[16], hc8K, 8, 16);
    if(*do_cluster_list[17] == 1) Write_Cluster_Xmol(f, *num_cluster_list[17], hc9A, 9, 17);
    if(*do_cluster_list[18] == 1) Write_Cluster_Xmol(f, *num_cluster_list[18], hc9B, 9, 18);
    if(*do_cluster_list[19] == 1) Write_Cluster_Xmol(f, *num_cluster_list[19], hc9K, 9, 19);
    if(*do_cluster_list[20] == 1) Write_Cluster_Xmol(f, *num_cluster_list[20], hc10A, 10, 20);
    if(*do_cluster_list[21] == 1) Write_Cluster_Xmol(f, *num_cluster_list[21], hc10B, 10, 21);
    if(*do_cluster_list[22] == 1) Write_Cluster_Xmol(f, *num_cluster_list[22], hc10K, 10, 22);
    if(*do_cluster_list[23] == 1) Write_Cluster_Xmol(f, *num_cluster_list[23], hc10W, 10, 23);
    if(*do_cluster_list[24] == 1) Write_Cluster_Xmol(f, *num_cluster_list[24], hc11A, 11, 24);
    if(*do_cluster_list[25] == 1) Write_Cluster_Xmol(f, *num_cluster_list[25], hc11B, 11, 25);
    if(*do_cluster_list[26] == 1) Write_Cluster_Xmol(f, *num_cluster_list[26], hc11C, 11, 26);
    if(*do_cluster_list[27] == 1) Write_Cluster_Xmol(f, *num_cluster_list[27], hc11E, 11, 27);
    if(*do_cluster_list[28] == 1) Write_Cluster_Xmol(f, *num_cluster_list[28], hc11F, 11, 28);
    if(*do_cluster_list[29] == 1) Write_Cluster_Xmol(f, *num_cluster_list[29], hc11W, 11, 29);
    if(*do_cluster_list[30] == 1) Write_Cluster_Xmol(f, *num_cluster_list[30], hc12A, 12, 30);
    if(*do_cluster_list[31] == 1) Write_Cluster_Xmol(f, *num_cluster_list[31], hc12B, 12, 31);
    if(*do_cluster_list[32] == 1) Write_Cluster_Xmol(f, *num_cluster_list[32], hc12D, 12, 32);
    if(*do_cluster_list[33] == 1) Write_Cluster_Xmol(f, *num_cluster_list[33], hc12E, 12, 33);
    if(*do_cluster_list[34] == 1) Write_Cluster_Xmol(f, *num_cluster_list[34], hc12K, 12, 34);
    if(*do_cluster_list[35] == 1) Write_Cluster_Xmol(f, *num_cluster_list[35], hc13A, 13, 35);
    if(*do_cluster_list[36] == 1) Write_Cluster_Xmol(f, *num_cluster_list[36], hc13B, 13, 36);
    if(*do_cluster_list[37] == 1) Write_Cluster_Xmol(f, *num_cluster_list[37], hc13K, 13, 37);
    if(*do_cluster_list[38] == 1) Write_Cluster_Xmol(f, *num_cluster_list[38], hcFCC, 13, 38);
    if(*do_cluster_list[39] == 1) Write_Cluster_Xmol(f, *num_cluster_list[39], hcHCP, 13, 39);
    if(*do_cluster_list[40] == 1) Write_Cluster_Xmol(f, *num_cluster_list[40], hcBCC_9, 9, 40);
    if(*do_cluster_list[41] == 1) Write_Cluster_Xmol(f, *num_cluster_list[41], hcBCC_15, 15, 41);
}

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
        fprintf(file_pointer,"%d",f);
        for(i=0; i<num_cluster_types; i++) {
            fprintf(file_pointer, "%.15lg	", pop_per_frame[i][f]);
        }
    }
    fclose(file_pointer);
}
