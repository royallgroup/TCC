#include "output.h"
#include "tools.h"
#include "globals.h"

void Write_Initialise_Raw_Files(int cluster_number);

void Write_initialise_Cluster_Files(int cluster_number);

void Write_Raw_Init() {
    int i;

    for(i=0; i<num_cluster_types; i++) {
        Write_Initialise_Raw_Files(i);
    }
}

void Write_Initialise_Raw_Files(int cluster_number) {
    char errMsg[1000];
    char output[1000];
    if(*do_cluster_list[cluster_number] == 1){
        sprintf(output, "%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_%s",
                fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
        raw_file_pointers[cluster_number] = fopen(output, "w");
        if (raw_file_pointers[cluster_number] == NULL) {
            sprintf(errMsg, "Write_Cluster_Init(): Error opening file %s", output);
            Error(errMsg);
        }
        fprintf(raw_file_pointers[cluster_number], "%s\n", output);
    }
    else {
        raw_file_pointers[cluster_number] = NULL;
    }
}

void Write_Raw_Xmol(int f, FILE *thefile, const char *sarr) {
    int i;

    fprintf(thefile,"%d\nframe %d of %d\n",N,f,TOTALFRAMES);
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

void Write_11A_cen_xmol(int f) {
    int i;

    int n11A_cen=0;

    for(i=0; i<N; i++) if(s11A_cen[i]=='O') ++n11A_cen;         // get total number of 13A centres

    fprintf(file_11A_cen_xmol,"%d\nframe %d of %d\n",n11A_cen,f,TOTALFRAMES);
    for(i=0; i<N; i++) {
        if (s11A_cen[i]=='O') fprintf(file_11A_cen_xmol ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
    }
}

void Write_13A_cen_xmol(int f) {
    int i;

    int n13A_cen=0;

    for(i=0; i<N; i++) if(s13A_cen[i]=='O') ++n13A_cen;         // get total number of 13A centres

    fprintf(file_13A_cen_xmol,"%d\nframe %d of %d\n",n13A_cen,f,FRAMES);
    for(i=0; i<N; i++) {
        if (s13A_cen[i]=='O') fprintf(file_13A_cen_xmol ,"O\t%.5lg\t%.5lg\t%.5lg\n", x[i], y[i], z[i]);
    }
}

void Write_Raw(int f) {
    int cluster_type;
    for(cluster_type=0; cluster_type<num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            Write_Raw_Xmol(f, raw_file_pointers[cluster_type], raw_cluster_list[cluster_type][0]);
        }
    }
}

void Write_Raw_Close() {
    int cluster_type;

    for(cluster_type=0; cluster_type< num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            fclose(raw_file_pointers[cluster_type]);
        }
    }
}

void Write_Cluster_Init() {
    int cluster_type;

    for(cluster_type=0; cluster_type<num_cluster_types; cluster_type++) {
        Write_initialise_Cluster_Files(cluster_type);
    }

}

void Write_initialise_Cluster_Files(int cluster_number) {
    char errMsg[1000];
    char output[1000];

    if(*do_cluster_list[cluster_number] == 1) {
        sprintf(output, "%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
                fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
        cluster_file_pointers[cluster_number] = fopen(output, "w");
        if (cluster_file_pointers[cluster_number] == NULL) {
            sprintf(errMsg, "Write_Cluster_Init(): Error opening file %s", output);
            Error(errMsg);
        }
        fprintf(cluster_file_pointers[cluster_number], "%s\n", output);
    }
    else {
        cluster_file_pointers[cluster_number] = NULL;
    }

}

void Write_Cluster_Close() {
    int cluster_type;

    for(cluster_type=0; cluster_type<num_cluster_types; cluster_type++) {
        if (*do_cluster_list[cluster_type] == 1) {
            fclose(cluster_file_pointers[cluster_type]);
        }
    }
}

void Write_Cluster_Xmol(int f, FILE *writeout, int n, int **hc, int clusSize) {
    int i,j;

    fprintf(writeout,"Frame Number%d\n",f);
    for (i=0;i<n;i++) {
        fprintf(writeout,"%d",hc[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hc[i][j]);
        fprintf(writeout,"	%d\n",hc[i][clusSize-1]);
    }
}

void Write_Cluster_sp3(int f, FILE *writeout) {
    int i,j;
    int clusSize=3;

    fprintf(writeout,"Frame Number%d\n",f);
    for (i=0;i<nsp3a;i++) {
        fprintf(writeout,"%d",hcsp3a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp3a[i][j]);
        fprintf(writeout,"	%d\n",hcsp3a[i][clusSize-1]);
    }
    for (i=0;i<nsp3b;i++) {
        fprintf(writeout,"%d",hcsp3b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp3b[i][j]);
        fprintf(writeout,"	%d\n",hcsp3b[i][clusSize-1]);
    }
    for (i=0;i<nsp3c;i++) {
        fprintf(writeout,"%d",hcsp3c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp3c[i][j]);
        fprintf(writeout,"	%d\n",hcsp3c[i][clusSize-1]);
    }
}

void Write_Cluster_sp4(int f, FILE *writeout) {
    int i,j;
    int clusSize=4;

    fprintf(writeout,"Frame Number%d\n",f);
    for (i=0;i<nsp4a;i++) {
        fprintf(writeout,"%d",hcsp4a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp4a[i][j]);
        fprintf(writeout,"	%d\n",hcsp4a[i][clusSize-1]);
    }
    for (i=0;i<nsp4b;i++) {
        fprintf(writeout,"%d",hcsp4b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp4b[i][j]);
        fprintf(writeout,"	%d\n",hcsp4b[i][clusSize-1]);
    }
    for (i=0;i<nsp4c;i++) {
        fprintf(writeout,"%d",hcsp4c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hcsp4c[i][j]);
        fprintf(writeout,"	%d\n",hcsp4c[i][clusSize-1]);
    }
}

void Write_Cluster_sp5(int f, FILE *writeout) {
    int i, j;
    int clusSize = 5;

    fprintf(writeout,"Frame Number%d\n",f);
    for (i = 0; i < nsp5a; i++) {
        fprintf(writeout, "%d", hcsp5a[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", hcsp5a[i][j]);
        fprintf(writeout, "	%d\n", hcsp5a[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5b; i++) {
        fprintf(writeout, "%d", hcsp5b[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", hcsp5b[i][j]);
        fprintf(writeout, "	%d\n", hcsp5b[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5c; i++) {
        fprintf(writeout, "%d", hcsp5c[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", hcsp5c[i][j]);
        fprintf(writeout, "	%d\n", hcsp5c[i][clusSize - 1]);
    }
}

void Write_Cluster(int f) {
    if(*do_cluster_list[0] == 1) Write_Cluster_sp3(f, raw_file_pointers[0]);
    if(*do_cluster_list[1] == 1) Write_Cluster_Xmol(f, raw_file_pointers[1], *num_cluster_list[1], hcsp3a, 3);
    if(*do_cluster_list[2] == 1) Write_Cluster_Xmol(f, raw_file_pointers[2], *num_cluster_list[2], hcsp3b, 4);
    if(*do_cluster_list[3] == 1) Write_Cluster_Xmol(f, raw_file_pointers[3], *num_cluster_list[3], hcsp3c, 5);
    if(*do_cluster_list[4] == 1) Write_Cluster_sp4(f, raw_file_pointers[4]);
    if(*do_cluster_list[5] == 1) Write_Cluster_Xmol(f, raw_file_pointers[5], *num_cluster_list[5], hcsp4a, 4);
    if(*do_cluster_list[6] == 1) Write_Cluster_Xmol(f, raw_file_pointers[6], *num_cluster_list[6], hcsp4b, 5);
    if(*do_cluster_list[7] == 1) Write_Cluster_Xmol(f, raw_file_pointers[7], *num_cluster_list[7], hcsp4c, 6);
    if(*do_cluster_list[8] == 1) Write_Cluster_sp5(f, raw_file_pointers[8]);
    if(*do_cluster_list[9] == 1) Write_Cluster_Xmol(f, raw_file_pointers[9], *num_cluster_list[9], hcsp5a, 5);
    if(*do_cluster_list[10] == 1) Write_Cluster_Xmol(f, raw_file_pointers[10], *num_cluster_list[10], hcsp5b, 6);
    if(*do_cluster_list[11] == 1) Write_Cluster_Xmol(f, raw_file_pointers[11], *num_cluster_list[11], hcsp5c, 7);
    if(*do_cluster_list[12] == 1) Write_Cluster_Xmol(f, raw_file_pointers[12], *num_cluster_list[12], hc6Z, 6);
    if(*do_cluster_list[13] == 1) Write_Cluster_Xmol(f, raw_file_pointers[13], *num_cluster_list[13], hc7K, 7);
    if(*do_cluster_list[14] == 1) Write_Cluster_Xmol(f, raw_file_pointers[14], *num_cluster_list[14], hc8A, 8);
    if(*do_cluster_list[15] == 1) Write_Cluster_Xmol(f, raw_file_pointers[15], *num_cluster_list[15], hc8B, 8);
    if(*do_cluster_list[16] == 1) Write_Cluster_Xmol(f, raw_file_pointers[16], *num_cluster_list[16], hc8K, 8);
    if(*do_cluster_list[17] == 1) Write_Cluster_Xmol(f, raw_file_pointers[17], *num_cluster_list[17], hc9A, 9);
    if(*do_cluster_list[18] == 1) Write_Cluster_Xmol(f, raw_file_pointers[18], *num_cluster_list[18], hc9B, 9);
    if(*do_cluster_list[19] == 1) Write_Cluster_Xmol(f, raw_file_pointers[19], *num_cluster_list[19], hc9K, 9);
    if(*do_cluster_list[20] == 1) Write_Cluster_Xmol(f, raw_file_pointers[20], *num_cluster_list[20], hc10A, 10);
    if(*do_cluster_list[21] == 1) Write_Cluster_Xmol(f, raw_file_pointers[21], *num_cluster_list[21], hc10B, 10);
    if(*do_cluster_list[22] == 1) Write_Cluster_Xmol(f, raw_file_pointers[22], *num_cluster_list[22], hc10K, 10);
    if(*do_cluster_list[23] == 1) Write_Cluster_Xmol(f, raw_file_pointers[23], *num_cluster_list[23], hc10W, 10);
    if(*do_cluster_list[24] == 1) Write_Cluster_Xmol(f, raw_file_pointers[24], *num_cluster_list[24], hc11A, 11);
    if(*do_cluster_list[25] == 1) Write_Cluster_Xmol(f, raw_file_pointers[25], *num_cluster_list[25], hc11B, 11);
    if(*do_cluster_list[26] == 1) Write_Cluster_Xmol(f, raw_file_pointers[26], *num_cluster_list[26], hc11C, 11);
    if(*do_cluster_list[27] == 1) Write_Cluster_Xmol(f, raw_file_pointers[27], *num_cluster_list[27], hc11E, 11);
    if(*do_cluster_list[28] == 1) Write_Cluster_Xmol(f, raw_file_pointers[28], *num_cluster_list[28], hc11F, 11);
    if(*do_cluster_list[29] == 1) Write_Cluster_Xmol(f, raw_file_pointers[29], *num_cluster_list[29], hc11W, 11);
    if(*do_cluster_list[30] == 1) Write_Cluster_Xmol(f, raw_file_pointers[30], *num_cluster_list[30], hc12A, 12);
    if(*do_cluster_list[31] == 1) Write_Cluster_Xmol(f, raw_file_pointers[31], *num_cluster_list[31], hc12B, 12);
    if(*do_cluster_list[32] == 1) Write_Cluster_Xmol(f, raw_file_pointers[32], *num_cluster_list[32], hc12D, 12);
    if(*do_cluster_list[33] == 1) Write_Cluster_Xmol(f, raw_file_pointers[33], *num_cluster_list[33], hc12E, 12);
    if(*do_cluster_list[34] == 1) Write_Cluster_Xmol(f, raw_file_pointers[34], *num_cluster_list[34], hc12K, 12);
    if(*do_cluster_list[35] == 1) Write_Cluster_Xmol(f, raw_file_pointers[35], *num_cluster_list[35], hc13A, 13);
    if(*do_cluster_list[36] == 1) Write_Cluster_Xmol(f, raw_file_pointers[36], *num_cluster_list[36], hc13B, 13);
    if(*do_cluster_list[37] == 1) Write_Cluster_Xmol(f, raw_file_pointers[37], *num_cluster_list[37], hc13K, 13);
    if(*do_cluster_list[38] == 1) Write_Cluster_Xmol(f, raw_file_pointers[38], *num_cluster_list[38], hcFCC, 13);
    if(*do_cluster_list[39] == 1) Write_Cluster_Xmol(f, raw_file_pointers[39], *num_cluster_list[39], hcHCP, 13);
    if(*do_cluster_list[40] == 1) Write_Cluster_Xmol(f, raw_file_pointers[40], *num_cluster_list[40], hcBCC_9, 9);
    if(*do_cluster_list[41] == 1) Write_Cluster_Xmol(f, raw_file_pointers[41], *num_cluster_list[41], hcBCC_15, 15);
}

void Write_Pop_Per_Frame(int f) {
    char errMsg[1000], output[1000];
    int i;

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    printf("Writing pop_per_frame %s\n",output);
    fPopPerFrame=fopen(output, "w");
    if (fPopPerFrame==NULL)  {
        sprintf(errMsg,"main() : Error opening file %s",output);	// Always test file open
        Error(errMsg);
    }
    fprintf(fPopPerFrame,"%s\n",output);

    fprintf(fPopPerFrame,"frame	");
    for(i=0; i<num_cluster_types; i++) {
        fprintf(fPopPerFrame, "%s	", cluster_names[i]);
    }
    fprintf(fPopPerFrame,"\n");

    fprintf(fPopPerFrame,"mean	");
    for(i=0; i<num_cluster_types; i++) {
        if(*do_cluster_list[i] == 1) {
            fprintf(fPopPerFrame, "%.15lg	", mean_pop_per_frame[i]);
        }
        else {
            fprintf(fPopPerFrame, "NA	");
        }
    }
    fprintf(fPopPerFrame,"\n");

    for (f=0;f<FRAMES;f++) {
        fprintf(fPopPerFrame,"%d",f);
        for(i=0; i<num_cluster_types; i++) {
            fprintf(fPopPerFrame, "%.15lg	", pop_per_frame[i][f]);
        }
    }
    fclose(fPopPerFrame);
    printf("Closed file %s\n\n",output);
}
