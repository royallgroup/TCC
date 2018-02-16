#include "output.h"
#include "tools.h"
#include "globals.h"

void Write_Initialise_Raw_Files(int cluster_number);

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

void Write_Raw_Xmol(int f, FILE *thefile, char *sarr) {
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
    if(*do_cluster_list[3] == 1) Write_Raw_Xmol(f,raw_file_pointers[3],&s5A[0]);
    if(*do_cluster_list[4] == 1) Write_Raw_Xmol(f,raw_file_pointers[4],&ssp4[0]);
    if(*do_cluster_list[5] == 1) Write_Raw_Xmol(f,raw_file_pointers[5],&ssp4a[0]);
    if(*do_cluster_list[6] == 1) Write_Raw_Xmol(f,raw_file_pointers[6],&ssp4b[0]);
    if(*do_cluster_list[7] == 1) Write_Raw_Xmol(f,raw_file_pointers[7],&s6A[0]);
    if(*do_cluster_list[8] == 1) Write_Raw_Xmol(f,raw_file_pointers[8],&s6Z[0]);
    if(*do_cluster_list[9] == 1) Write_Raw_Xmol(f,raw_file_pointers[9],&s7K[0]);
    if(*do_cluster_list[10] == 1) Write_Raw_Xmol(f,raw_file_pointers[10],&ssp5[0]);
    if(*do_cluster_list[11] == 1) Write_Raw_Xmol(f,raw_file_pointers[11],&ssp5a[0]);
    if(*do_cluster_list[12] == 1) Write_Raw_Xmol(f,raw_file_pointers[12],&ssp5b[0]);
    if(*do_cluster_list[13] == 1) Write_Raw_Xmol(f,raw_file_pointers[13],&s7A[0]);
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
    char errMsg[1000];
    char output[1000];

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp3=fopen(output,"w");
    if (wsp3==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp3,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp3a=fopen(output,"w");
    if (wsp3a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp3a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp3b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp3b=fopen(output,"w");
    if (wsp3b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp3b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_5A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w5A=fopen(output,"w");
    if (w5A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w5A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp4=fopen(output,"w");
    if (wsp4==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp4,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp4a=fopen(output,"w");
    if (wsp4a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp4a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp4b=fopen(output,"w");
    if (wsp4b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp4b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp4c",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp4c=fopen(output,"w");
    if (wsp4c==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp4c,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_6A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w6A=fopen(output,"w");
    if (w6A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w6A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_6Z",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w6Z=fopen(output,"w");
    if (w6Z==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w6Z,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp5=fopen(output,"w");
    if (wsp5==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp5,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp5a=fopen(output,"w");
    if (wsp5a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp5a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_sp5b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wsp5b=fopen(output,"w");
    if (wsp5b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wsp5b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_7A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w7A=fopen(output,"w");
    if (w7A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w7A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_7K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w7K=fopen(output,"w");
    if (w7K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w7K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w8A=fopen(output,"w");
    if (w8A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w8A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w8B=fopen(output,"w");
    if (w8B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w8B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_8K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w8K=fopen(output,"w");
    if (w8K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w8K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w9A=fopen(output,"w");
    if (w9A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w9A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w9B=fopen(output,"w");
    if (w9B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w9B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_9K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w9K=fopen(output,"w");
    if (w9K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w9K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w10A=fopen(output,"w");
    if (w10A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w10A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w10B=fopen(output,"w");
    if (w10B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w10B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w10K=fopen(output,"w");
    if (w10K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w10K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_10W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w10W=fopen(output,"w");
    if (w10W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w10W,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11A=fopen(output,"w");
    if (w11A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11B=fopen(output,"w");
    if (w11B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11C",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11C=fopen(output,"w");
    if (w11C==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11C,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11E=fopen(output,"w");
    if (w11E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11E,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11F",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11F=fopen(output,"w");
    if (w11F==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11F,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_11W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w11W=fopen(output,"w");
    if (w11W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w11W,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w12A=fopen(output,"w");
    if (w12A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w12A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w12B=fopen(output,"w");
    if (w12B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w12B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12D",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w12D=fopen(output,"w");
    if (w12D==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w12D,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w12E=fopen(output,"w");
    if (w12E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w12E,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_12K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w12K=fopen(output,"w");
    if (w12K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w12K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w13A=fopen(output,"w");
    if (w13A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w13A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w13B=fopen(output,"w");
    if (w13B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w13B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_13K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    w13K=fopen(output,"w");
    if (w13K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(w13K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_FCC",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wFCC=fopen(output,"w");
    if (wFCC==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wFCC,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_HCP",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wHCP=fopen(output,"w");
    if (wHCP==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wHCP,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_BCC_9",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wBCC_9=fopen(output,"w");
    if (wBCC_9==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wBCC_9,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_BCC_15",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    wBCC_15=fopen(output,"w");
    if (wBCC_15==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(wBCC_15,"%s\n",output);
}

void Write_Cluster_Close() {
    fclose(wsp3);
    fclose(wsp3a);
    fclose(wsp3b);
    fclose(w5A);
    fclose(wsp4);
    fclose(wsp4a);
    fclose(wsp4b);
    fclose(wsp4c);
    fclose(w6A);
    fclose(w6Z);
    fclose(w7K);
    fclose(wsp5);
    fclose(wsp5a);
    fclose(wsp5b);
    fclose(w7A);
    fclose(w8A);
    fclose(w8B);
    fclose(w8K);
    fclose(w9A);
    fclose(w9B);
    fclose(w9K);
    fclose(w10A);
    fclose(w10B);
    fclose(w10K);
    fclose(w10W);
    fclose(w11A);
    fclose(w11B);
    fclose(w11C);
    fclose(w11E);
    fclose(w11F);
    fclose(w11W);
    fclose(w12A);
    fclose(w12B);
    fclose(w12D);
    fclose(w12E);
    fclose(w12K);
    fclose(w13A);
    fclose(w13B);
    fclose(w13K);
    fclose(wFCC);
    fclose(wHCP);
    fclose(wBCC_9);
    fclose(wBCC_15);
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
    Write_Cluster_sp3(f, wsp3);
    Write_Cluster_Xmol(f, wsp3a, nsp3a, sp3a, 3);
    Write_Cluster_Xmol(f, wsp3b, nsp3b, sp3b, 4);
    Write_Cluster_Xmol(f, w5A, nsp3c, sp3c, 5);
    Write_Cluster_sp4(f, wsp4);
    Write_Cluster_Xmol(f, wsp4a, nsp4a, sp4a, 4);
    Write_Cluster_Xmol(f, wsp4b, nsp4b, sp4b, 5);
    Write_Cluster_Xmol(f, wsp4c, nsp4c, sp4c, 6);
    Write_Cluster_Xmol(f, w6A, n6A, hc6A, 6);
    Write_Cluster_sp5(f, wsp5);
    Write_Cluster_Xmol(f, wsp5a, nsp5a, sp5a, 5);
    Write_Cluster_Xmol(f, wsp5b, nsp5b, sp5b, 6);
    Write_Cluster_Xmol(f, w7A, nsp5c, sp5c, 7);
    Write_Cluster_Xmol(f, w6Z, n6Z, hc6Z, 6);
    Write_Cluster_Xmol(f, w7K, n7K, hc7K, 7);
    Write_Cluster_Xmol(f, w8A, n8A, hc8A, 8);
    Write_Cluster_Xmol(f, w8B, n8B, hc8B, 8);
    Write_Cluster_Xmol(f, w8K, n8K, hc8K, 8);
    Write_Cluster_Xmol(f, w9A, n9A, hc9A, 9);
    Write_Cluster_Xmol(f, w9B, n9B, hc9B, 9);
    Write_Cluster_Xmol(f, w9K, n9K, hc9K, 9);
    Write_Cluster_Xmol(f, w10A, n10A, hc10A, 10);
    Write_Cluster_Xmol(f, w10B, n10B, hc10B, 10);
    Write_Cluster_Xmol(f, w10K, n10K, hc10K, 10);
    Write_Cluster_Xmol(f, w10W, n10W, hc10W, 10);
    Write_Cluster_Xmol(f, w11A, n11A, hc11A, 11);
    Write_Cluster_Xmol(f, w11B, n11B, hc11B, 11);
    Write_Cluster_Xmol(f, w11C, n11C, hc11C, 11);
    Write_Cluster_Xmol(f, w11E, n11E, hc11E, 11);
    Write_Cluster_Xmol(f, w11F, n11F, hc11F, 11);
    Write_Cluster_Xmol(f, w11W, n11W, hc11W, 11);
    Write_Cluster_Xmol(f, w12A, n12A, hc12A, 12);
    Write_Cluster_Xmol(f, w12B, n12B, hc12B, 12);
    Write_Cluster_Xmol(f, w12D, n12D, hc12D, 12);
    Write_Cluster_Xmol(f, w12E, n12E, hc12E, 12);
    Write_Cluster_Xmol(f, w12K, n12K, hc12K, 12);
    Write_Cluster_Xmol(f, w13A, n13A, hc13A, 13);
    Write_Cluster_Xmol(f, w13B, n13B, hc13B, 13);
    Write_Cluster_Xmol(f, w13K, n13K, hc13K, 13);
    Write_Cluster_Xmol(f, wFCC, nFCC, hcFCC, 13);
    Write_Cluster_Xmol(f, wHCP, nHCP, hcHCP, 13);
    Write_Cluster_Xmol(f, wBCC_9, nBCC_9, hcBCC_9, 9);
    Write_Cluster_Xmol(f, wBCC_15, nBCC_15, hcBCC_15, 15);
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


    fprintf(fPopPerFrame,"mean	-	-	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",mean_pop_per_frame_sp3,mean_pop_per_frame_sp3a,mean_pop_per_frame_sp3b,mean_pop_per_frame_sp3c,mean_pop_per_frame_sp4,mean_pop_per_frame_sp4a,mean_pop_per_frame_sp4b,mean_pop_per_frame_6A,mean_pop_per_frame_6Z,mean_pop_per_frame_sp5,mean_pop_per_frame_sp5a,mean_pop_per_frame_sp5b,mean_pop_per_frame_sp5c,mean_pop_per_frame_7K,mean_pop_per_frame_8A,mean_pop_per_frame_8B,mean_pop_per_frame_8K,mean_pop_per_frame_9A,mean_pop_per_frame_9B,mean_pop_per_frame_9K,mean_pop_per_frame_10A,mean_pop_per_frame_10B,mean_pop_per_frame_10K,mean_pop_per_frame_10W);
    fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",mean_pop_per_frame_11A,mean_pop_per_frame_11B,mean_pop_per_frame_11C,mean_pop_per_frame_11E,mean_pop_per_frame_11F,mean_pop_per_frame_11W,mean_pop_per_frame_12A,mean_pop_per_frame_12B,mean_pop_per_frame_12D,mean_pop_per_frame_12E,mean_pop_per_frame_12K,mean_pop_per_frame_13A,mean_pop_per_frame_13B,mean_pop_per_frame_13K,mean_pop_per_frame_FCC,mean_pop_per_frame_HCP,mean_pop_per_frame_BCC_9,mean_pop_per_frame_BCC_15);
    for (f=0;f<FRAMES;f++) {
        fprintf(fPopPerFrame,"%d	%.15lg",f,(double)f*FRAMETSTEP*SAMPLEFREQ);
        fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg",pop_per_frame_sp3[f],pop_per_frame_sp3a[f],pop_per_frame_sp3b[f],pop_per_frame_sp3c[f],pop_per_frame_sp4[f],pop_per_frame_sp4a[f],pop_per_frame_sp4b[f],pop_per_frame_6A[f],pop_per_frame_6Z[f],pop_per_frame_sp5[f],pop_per_frame_sp5a[f],pop_per_frame_sp5b[f],pop_per_frame_sp5c[f],pop_per_frame_7K[f],pop_per_frame_8A[f],pop_per_frame_8B[f],pop_per_frame_8K[f],pop_per_frame_9A[f],pop_per_frame_9B[f],pop_per_frame_9K[f],pop_per_frame_10A[f],pop_per_frame_10B[f],pop_per_frame_10K[f],pop_per_frame_10W[f]);
        fprintf(fPopPerFrame,"	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg	%.15lg\n",pop_per_frame_11A[f],pop_per_frame_11B[f],pop_per_frame_11C[f],pop_per_frame_11E[f],pop_per_frame_11F[f],pop_per_frame_11W[f],pop_per_frame_12A[f],pop_per_frame_12B[f],pop_per_frame_12D[f],pop_per_frame_12E[f],pop_per_frame_12K[f],pop_per_frame_13A[f],pop_per_frame_13B[f],pop_per_frame_13K[f],pop_per_frame_FCC[f],pop_per_frame_HCP[f],pop_per_frame_BCC_9[f],pop_per_frame_BCC_15[f]);
    }
    fclose(fPopPerFrame);
    printf("Closed file %s\n\n",output);
}
