#include "output.h"
#include "tools.h"
#include "globals.h"

void Write_Initialise_Raw_Files(int cluster_number);

void Write_initialise_Cluster_Files(int cluster_number);

void Write_Raw_Init() {
    int i;

    raw_file_pointers = malloc(sizeof(FILE*) * num_cluster_types);

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
    if(*do_cluster_list[0] == 1) Write_Raw_Xmol(f,raw_file_pointers[0],&ssp3[0]);
    if(*do_cluster_list[1] == 1) Write_Raw_Xmol(f,raw_file_pointers[1],&ssp3a[0]);
    if(*do_cluster_list[2] == 1) Write_Raw_Xmol(f,raw_file_pointers[2],&ssp3b[0]);
    if(*do_cluster_list[3] == 1) Write_Raw_Xmol(f,raw_file_pointers[3],&ssp3c[0]);
    if(*do_cluster_list[4] == 1) Write_Raw_Xmol(f,raw_file_pointers[4],&ssp4[0]);
    if(*do_cluster_list[5] == 1) Write_Raw_Xmol(f,raw_file_pointers[5],&ssp4a[0]);
    if(*do_cluster_list[6] == 1) Write_Raw_Xmol(f,raw_file_pointers[6],&ssp4b[0]);
    if(*do_cluster_list[7] == 1) Write_Raw_Xmol(f,raw_file_pointers[7],&ssp4c[0]);
    if(*do_cluster_list[8] == 1) Write_Raw_Xmol(f,raw_file_pointers[8],&s6Z[0]);
    if(*do_cluster_list[9] == 1) Write_Raw_Xmol(f,raw_file_pointers[9],&s7K[0]);
    if(*do_cluster_list[10] == 1) Write_Raw_Xmol(f,raw_file_pointers[10],&ssp5[0]);
    if(*do_cluster_list[11] == 1) Write_Raw_Xmol(f,raw_file_pointers[11],&ssp5a[0]);
    if(*do_cluster_list[12] == 1) Write_Raw_Xmol(f,raw_file_pointers[12],&ssp5b[0]);
    if(*do_cluster_list[13] == 1) Write_Raw_Xmol(f,raw_file_pointers[13],&ssp5c[0]);
    if(*do_cluster_list[14] == 1) Write_Raw_Xmol(f,raw_file_pointers[14],&s8A[0]);
    if(*do_cluster_list[15] == 1) Write_Raw_Xmol(f,raw_file_pointers[15],&s8B[0]);
    if(*do_cluster_list[16] == 1) Write_Raw_Xmol(f,raw_file_pointers[16],&s8K[0]);
    if(*do_cluster_list[17] == 1) Write_Raw_Xmol(f,raw_file_pointers[17],&s9A[0]);
    if(*do_cluster_list[18] == 1) Write_Raw_Xmol(f,raw_file_pointers[18],&s9B[0]);
    if(*do_cluster_list[19] == 1) Write_Raw_Xmol(f,raw_file_pointers[19],&s9K[0]);
    if(*do_cluster_list[20] == 1) Write_Raw_Xmol(f,raw_file_pointers[20],&s10A[0]);
    if(*do_cluster_list[21] == 1) Write_Raw_Xmol(f,raw_file_pointers[21],&s10B[0]);
    if(*do_cluster_list[22] == 1) Write_Raw_Xmol(f,raw_file_pointers[22],&s10K[0]);
    if(*do_cluster_list[23] == 1) Write_Raw_Xmol(f,raw_file_pointers[23],&s10W[0]);
    if(*do_cluster_list[24] == 1) Write_Raw_Xmol(f,raw_file_pointers[24],&s11A[0]);
    if(*do_cluster_list[25] == 1) Write_Raw_Xmol(f,raw_file_pointers[25],&s11B[0]);
    if(*do_cluster_list[26] == 1) Write_Raw_Xmol(f,raw_file_pointers[26],&s11C[0]);
    if(*do_cluster_list[27] == 1) Write_Raw_Xmol(f,raw_file_pointers[27],&s11E[0]);
    if(*do_cluster_list[28] == 1) Write_Raw_Xmol(f,raw_file_pointers[28],&s11F[0]);
    if(*do_cluster_list[29] == 1) Write_Raw_Xmol(f,raw_file_pointers[29],&s11W[0]);
    if(*do_cluster_list[30] == 1) Write_Raw_Xmol(f,raw_file_pointers[30],&s12A[0]);
    if(*do_cluster_list[31] == 1) Write_Raw_Xmol(f,raw_file_pointers[31],&s12B[0]);
    if(*do_cluster_list[32] == 1) Write_Raw_Xmol(f,raw_file_pointers[32],&s12D[0]);
    if(*do_cluster_list[33] == 1) Write_Raw_Xmol(f,raw_file_pointers[33],&s12E[0]);
    if(*do_cluster_list[34] == 1) Write_Raw_Xmol(f,raw_file_pointers[34],&s12K[0]);
    if(*do_cluster_list[35] == 1) Write_Raw_Xmol(f,raw_file_pointers[35],&s13A[0]);
    if(*do_cluster_list[36] == 1) Write_Raw_Xmol(f,raw_file_pointers[36],&s13B[0]);
    if(*do_cluster_list[37] == 1) Write_Raw_Xmol(f,raw_file_pointers[37],&s13K[0]);
    if(*do_cluster_list[38] == 1) Write_Raw_Xmol(f,raw_file_pointers[38],&sFCC[0]);
    if(*do_cluster_list[39] == 1) Write_Raw_Xmol(f,raw_file_pointers[39],&sHCP[0]);
    if(*do_cluster_list[40] == 1) Write_Raw_Xmol(f,raw_file_pointers[40],&sBCC_9[0]);
    if(*do_cluster_list[41] == 1) Write_Raw_Xmol(f,raw_file_pointers[41],&sBCC_15[0]);
}

void Write_Raw_Close() {
    int i;

    for(i=0; i< num_cluster_types; i++) {
        fclose(raw_file_pointers[i]);
    }
}

void Write_Cluster_Init() {
    int i;

    cluster_file_pointers = malloc(sizeof(FILE*) * num_cluster_types);

    for(i=0; i<num_cluster_types; i++) {
        Write_initialise_Cluster_Files(i);
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
    int i;

    for(i=0; i<num_cluster_types; i++) {
        fclose(cluster_file_pointers[i]);
    }
}

void Write_Cluster_Xmol(int f, FILE *writeout, int *n, int **hc, int clusSize) {
    int i,j;

    fprintf(writeout,"%d\n",n[f]);
    for (i=0;i<n[f];i++) {
        fprintf(writeout,"%d",hc[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",hc[i][j]);
        fprintf(writeout,"	%d\n",hc[i][clusSize-1]);
    }
}

void Write_Cluster_sp3(int f, FILE *writeout) {
    int i,j;
    int clusSize=3;

    fprintf(writeout,"%d\n",nsp3[f]);
    for (i=0;i<nsp3a[f];i++) {
        fprintf(writeout,"%d",sp3a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3a[i][j]);
        fprintf(writeout,"	%d\n",sp3a[i][clusSize-1]);
    }
    for (i=0;i<nsp3b[f];i++) {
        fprintf(writeout,"%d",sp3b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3b[i][j]);
        fprintf(writeout,"	%d\n",sp3b[i][clusSize-1]);
    }
    for (i=0;i<nsp3c[f];i++) {
        fprintf(writeout,"%d",sp3c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp3c[i][j]);
        fprintf(writeout,"	%d\n",sp3c[i][clusSize-1]);
    }
}

void Write_Cluster_sp4(int f, FILE *writeout) {
    int i,j;
    int clusSize=4;

    fprintf(writeout,"%d\n",nsp4[f]);
    for (i=0;i<nsp4a[f];i++) {
        fprintf(writeout,"%d",sp4a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4a[i][j]);
        fprintf(writeout,"	%d\n",sp4a[i][clusSize-1]);
    }
    for (i=0;i<nsp4b[f];i++) {
        fprintf(writeout,"%d",sp4b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4b[i][j]);
        fprintf(writeout,"	%d\n",sp4b[i][clusSize-1]);
    }
    for (i=0;i<nsp4c[f];i++) {
        fprintf(writeout,"%d",sp4c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp4c[i][j]);
        fprintf(writeout,"	%d\n",sp4c[i][clusSize-1]);
    }
}

void Write_Cluster_sp5(int f, FILE *writeout) {
    int i, j;
    int clusSize = 5;

    fprintf(writeout, "%d\n", nsp5[f]);
    for (i = 0; i < nsp5a[f]; i++) {
        fprintf(writeout, "%d", sp5a[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", sp5a[i][j]);
        fprintf(writeout, "	%d\n", sp5a[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5b[f]; i++) {
        fprintf(writeout, "%d", sp5b[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", sp5b[i][j]);
        fprintf(writeout, "	%d\n", sp5b[i][clusSize - 1]);
    }
    for (i = 0; i < nsp5c[f]; i++) {
        fprintf(writeout, "%d", sp5c[i][0]);
        for (j = 1; j < clusSize - 1; j++) fprintf(writeout, "	%d", sp5c[i][j]);
        fprintf(writeout, "	%d\n", sp5c[i][clusSize - 1]);
    }
}

void Write_Cluster(int f) {
    if(*do_cluster_list[0] == 1) Write_Cluster_sp3(f, raw_file_pointers[0]);
    if(*do_cluster_list[1] == 1) Write_Cluster_Xmol(f, raw_file_pointers[1], nsp3a, sp3a, 3);
    if(*do_cluster_list[2] == 1) Write_Cluster_Xmol(f, raw_file_pointers[2], nsp3b, sp3b, 4);
    if(*do_cluster_list[3] == 1) Write_Cluster_Xmol(f, raw_file_pointers[3], nsp3c, sp3c, 5);
    if(*do_cluster_list[4] == 1) Write_Cluster_sp4(f, raw_file_pointers[4]);
    if(*do_cluster_list[5] == 1) Write_Cluster_Xmol(f, raw_file_pointers[5], nsp4a, sp4a, 4);
    if(*do_cluster_list[6] == 1) Write_Cluster_Xmol(f, raw_file_pointers[6], nsp4b, sp4b, 5);
    if(*do_cluster_list[7] == 1) Write_Cluster_Xmol(f, raw_file_pointers[7], nsp4c, sp4c, 6);
    if(*do_cluster_list[8] == 1) Write_Cluster_sp5(f, raw_file_pointers[8]);
    if(*do_cluster_list[9] == 1) Write_Cluster_Xmol(f, raw_file_pointers[9], nsp5a, sp5a, 5);
    if(*do_cluster_list[10] == 1) Write_Cluster_Xmol(f, raw_file_pointers[10], nsp5b, sp5b, 6);
    if(*do_cluster_list[11] == 1) Write_Cluster_Xmol(f, raw_file_pointers[11], nsp5c, sp5c, 7);
    if(*do_cluster_list[12] == 1) Write_Cluster_Xmol(f, raw_file_pointers[12], n6Z, hc6Z, 6);
    if(*do_cluster_list[13] == 1) Write_Cluster_Xmol(f, raw_file_pointers[13], n7K, hc7K, 7);
    if(*do_cluster_list[14] == 1) Write_Cluster_Xmol(f, raw_file_pointers[14], n8A, hc8A, 8);
    if(*do_cluster_list[15] == 1) Write_Cluster_Xmol(f, raw_file_pointers[15], n8B, hc8B, 8);
    if(*do_cluster_list[16] == 1) Write_Cluster_Xmol(f, raw_file_pointers[16], n8K, hc8K, 8);
    if(*do_cluster_list[17] == 1) Write_Cluster_Xmol(f, raw_file_pointers[17], n9A, hc9A, 9);
    if(*do_cluster_list[18] == 1) Write_Cluster_Xmol(f, raw_file_pointers[18], n9B, hc9B, 9);
    if(*do_cluster_list[19] == 1) Write_Cluster_Xmol(f, raw_file_pointers[19], n9K, hc9K, 9);
    if(*do_cluster_list[20] == 1) Write_Cluster_Xmol(f, raw_file_pointers[20], n10A, hc10A, 10);
    if(*do_cluster_list[21] == 1) Write_Cluster_Xmol(f, raw_file_pointers[21], n10B, hc10B, 10);
    if(*do_cluster_list[22] == 1) Write_Cluster_Xmol(f, raw_file_pointers[22], n10K, hc10K, 10);
    if(*do_cluster_list[23] == 1) Write_Cluster_Xmol(f, raw_file_pointers[23], n10W, hc10W, 10);
    if(*do_cluster_list[24] == 1) Write_Cluster_Xmol(f, raw_file_pointers[24], n11A, hc11A, 11);
    if(*do_cluster_list[25] == 1) Write_Cluster_Xmol(f, raw_file_pointers[25], n11B, hc11B, 11);
    if(*do_cluster_list[26] == 1) Write_Cluster_Xmol(f, raw_file_pointers[26], n11C, hc11C, 11);
    if(*do_cluster_list[27] == 1) Write_Cluster_Xmol(f, raw_file_pointers[27], n11E, hc11E, 11);
    if(*do_cluster_list[28] == 1) Write_Cluster_Xmol(f, raw_file_pointers[28], n11F, hc11F, 11);
    if(*do_cluster_list[29] == 1) Write_Cluster_Xmol(f, raw_file_pointers[29], n11W, hc11W, 11);
    if(*do_cluster_list[30] == 1) Write_Cluster_Xmol(f, raw_file_pointers[30], n12A, hc12A, 12);
    if(*do_cluster_list[31] == 1) Write_Cluster_Xmol(f, raw_file_pointers[31], n12B, hc12B, 12);
    if(*do_cluster_list[32] == 1) Write_Cluster_Xmol(f, raw_file_pointers[32], n12D, hc12D, 12);
    if(*do_cluster_list[33] == 1) Write_Cluster_Xmol(f, raw_file_pointers[33], n12E, hc12E, 12);
    if(*do_cluster_list[34] == 1) Write_Cluster_Xmol(f, raw_file_pointers[34], n12K, hc12K, 12);
    if(*do_cluster_list[35] == 1) Write_Cluster_Xmol(f, raw_file_pointers[35], n13A, hc13A, 13);
    if(*do_cluster_list[36] == 1) Write_Cluster_Xmol(f, raw_file_pointers[36], n13B, hc13B, 13);
    if(*do_cluster_list[37] == 1) Write_Cluster_Xmol(f, raw_file_pointers[37], n13K, hc13K, 13);
    if(*do_cluster_list[38] == 1) Write_Cluster_Xmol(f, raw_file_pointers[38], nFCC, hcFCC, 13);
    if(*do_cluster_list[39] == 1) Write_Cluster_Xmol(f, raw_file_pointers[39], nHCP, hcHCP, 13);
    if(*do_cluster_list[40] == 1) Write_Cluster_Xmol(f, raw_file_pointers[40], nBCC_9, hcBCC_9, 9);
    if(*do_cluster_list[41] == 1) Write_Cluster_Xmol(f, raw_file_pointers[41], nBCC_15, hcBCC_15, 15);
}

void Write_Pop_Per_Frame(int f) {
    char errMsg[1000], output[1000];

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    printf("Writing pop_per_frame %s\n",output);
    fPopPerFrame=fopen(output, "w");
    if (fPopPerFrame==NULL)  {
        sprintf(errMsg,"main() : Error opening file %s",output);	// Always test file open
        Error(errMsg);
    }
    fprintf(fPopPerFrame,"%s\n",output);

    fprintf(fPopPerFrame,"frame	time	time_norm_t_a	sp3	sp3a	sp3b	5A	sp4	sp4a	sp4b	6A	6Z	sp5	sp5a	sp5b	7A	7K	8A	8B	8K	9A	9B	9K	10A	10B	10K	10W");
    fprintf(fPopPerFrame,"	11A	11B	11C	11E	11F	11W	12A	12B	12D	12E	12K	13A	13B	13K	FCC	HCP	BCC_9	BCC_15\n");


    fprintf(fPopPerFrame,"mean	-	-	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",mean_pop_per_frame_sp3,mean_pop_per_frame_sp3a,mean_pop_per_frame_sp3b,mean_pop_per_frame_sp3c,mean_pop_per_frame_sp4,mean_pop_per_frame_sp4a,mean_pop_per_frame_sp4b,mean_pop_per_frame_sp4c,mean_pop_per_frame_6Z,mean_pop_per_frame_sp5,mean_pop_per_frame_sp5a,mean_pop_per_frame_sp5b,mean_pop_per_frame_sp5c,mean_pop_per_frame_7K,mean_pop_per_frame_8A,mean_pop_per_frame_8B,mean_pop_per_frame_8K,mean_pop_per_frame_9A,mean_pop_per_frame_9B,mean_pop_per_frame_9K,mean_pop_per_frame_10A,mean_pop_per_frame_10B,mean_pop_per_frame_10K,mean_pop_per_frame_10W);
    fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",mean_pop_per_frame_11A,mean_pop_per_frame_11B,mean_pop_per_frame_11C,mean_pop_per_frame_11E,mean_pop_per_frame_11F,mean_pop_per_frame_11W,mean_pop_per_frame_12A,mean_pop_per_frame_12B,mean_pop_per_frame_12D,mean_pop_per_frame_12E,mean_pop_per_frame_12K,mean_pop_per_frame_13A,mean_pop_per_frame_13B,mean_pop_per_frame_13K,mean_pop_per_frame_FCC,mean_pop_per_frame_HCP,mean_pop_per_frame_BCC_9,mean_pop_per_frame_BCC_15);
    for (f=0;f<FRAMES;f++) {
        fprintf(fPopPerFrame,"%d	%.15lg",f,(double)f*FRAMETSTEP*SAMPLEFREQ);
        fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",pop_per_frame_sp3[f],pop_per_frame_sp3a[f],pop_per_frame_sp3b[f],pop_per_frame_sp3c[f],pop_per_frame_sp4[f],pop_per_frame_sp4a[f],pop_per_frame_sp4b[f],pop_per_frame_sp4c[f],pop_per_frame_6Z[f],pop_per_frame_sp5[f],pop_per_frame_sp5a[f],pop_per_frame_sp5b[f],pop_per_frame_sp5c[f],pop_per_frame_7K[f],pop_per_frame_8A[f],pop_per_frame_8B[f],pop_per_frame_8K[f],pop_per_frame_9A[f],pop_per_frame_9B[f],pop_per_frame_9K[f],pop_per_frame_10A[f],pop_per_frame_10B[f],pop_per_frame_10K[f],pop_per_frame_10W[f]);
        fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",pop_per_frame_11A[f],pop_per_frame_11B[f],pop_per_frame_11C[f],pop_per_frame_11E[f],pop_per_frame_11F[f],pop_per_frame_11W[f],pop_per_frame_12A[f],pop_per_frame_12B[f],pop_per_frame_12D[f],pop_per_frame_12E[f],pop_per_frame_12K[f],pop_per_frame_13A[f],pop_per_frame_13B[f],pop_per_frame_13K[f],pop_per_frame_FCC[f],pop_per_frame_HCP[f],pop_per_frame_BCC_9[f],pop_per_frame_BCC_15[f]);
    }
    fclose(fPopPerFrame);
    printf("Closed file %s\n\n",output);
}
