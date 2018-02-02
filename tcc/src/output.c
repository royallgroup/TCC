#include "../inc/output.h"
#include "globals.h"
#include "tools.h"


void Write_Raw_Init() {
    char errMsg[1000];
    char output[1000];

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp3=fopen(output,"w");
    if (file_raw_sp3==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp3,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp3a=fopen(output,"w");
    if (file_raw_sp3a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp3a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp3b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp3b=fopen(output,"w");
    if (file_raw_sp3b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp3b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_5A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_5A=fopen(output,"w");
    if (file_raw_5A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_5A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp4=fopen(output,"w");
    if (file_raw_sp4==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp4,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp4a=fopen(output,"w");
    if (file_raw_sp4a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp4a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp4b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp4b=fopen(output,"w");
    if (file_raw_sp4b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp4b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_6A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_6A=fopen(output,"w");
    if (file_raw_6A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_6A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_6Z",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_6Z=fopen(output,"w");
    if (file_raw_6Z==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_6Z,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp5=fopen(output,"w");
    if (file_raw_sp5==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp5,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5a",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp5a=fopen(output,"w");
    if (file_raw_sp5a==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp5a,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_sp5b",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_sp5b=fopen(output,"w");
    if (file_raw_sp5b==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_sp5b,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_7A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_7A=fopen(output,"w");
    if (file_raw_7A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_7A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_7K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_7K=fopen(output,"w");
    if (file_raw_7K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_7K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_8A=fopen(output,"w");
    if (file_raw_8A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_8A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_8B=fopen(output,"w");
    if (file_raw_8B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_8B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_8K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_8K=fopen(output,"w");
    if (file_raw_8K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_8K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_9A=fopen(output,"w");
    if (file_raw_9A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_9A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_9B=fopen(output,"w");
    if (file_raw_9B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_9B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_9K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_9K=fopen(output,"w");
    if (file_raw_9K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_9K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_10A=fopen(output,"w");
    if (file_raw_10A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_10A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_10B=fopen(output,"w");
    if (file_raw_10B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_10B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_10K=fopen(output,"w");
    if (file_raw_10K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_10K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_10W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_10W=fopen(output,"w");
    if (file_raw_10W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_10W,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11A=fopen(output,"w");
    if (file_raw_11A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11B=fopen(output,"w");
    if (file_raw_11B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11C",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11C=fopen(output,"w");
    if (file_raw_11C==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11C,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11E=fopen(output,"w");
    if (file_raw_11E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11E,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11F",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11F=fopen(output,"w");
    if (file_raw_11F==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11F,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11W",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_11W=fopen(output,"w");
    if (file_raw_11W==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_11W,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_12A=fopen(output,"w");
    if (file_raw_12A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_12A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_12B=fopen(output,"w");
    if (file_raw_12B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_12B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12D",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_12D=fopen(output,"w");
    if (file_raw_12D==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_12D,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12E",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_12E=fopen(output,"w");
    if (file_raw_12E==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_12E,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_12K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_12K=fopen(output,"w");
    if (file_raw_12K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_12K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13A",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_13A=fopen(output,"w");
    if (file_raw_13A==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_13A,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13B",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_13B=fopen(output,"w");
    if (file_raw_13B==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_13B,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13K",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_13K=fopen(output,"w");
    if (file_raw_13K==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_13K,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_FCC",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_FCC=fopen(output,"w");
    if (file_raw_FCC==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_FCC,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_HCP",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_HCP=fopen(output,"w");
    if (file_raw_HCP==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_HCP,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_BCC_9",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_BCC_9=fopen(output,"w");
    if (file_raw_BCC_9==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_BCC_9,"%s\n",output);

    sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_BCC_15",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
    file_raw_BCC_15=fopen(output,"w");
    if (file_raw_BCC_15==NULL) { sprintf(errMsg,"Write_Cluster_Init(): Error opening file %s",output); Error(errMsg); }
    fprintf(file_raw_BCC_15,"%s\n",output);
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
    Write_Raw_Xmol(f,file_raw_sp3,&ssp3[0]);
    Write_Raw_Xmol(f,file_raw_sp3a,&ssp3a[0]);
    Write_Raw_Xmol(f,file_raw_sp3b,&ssp3b[0]);
    Write_Raw_Xmol(f,file_raw_5A,&s5A[0]);
    Write_Raw_Xmol(f,file_raw_sp4,&ssp4[0]);
    Write_Raw_Xmol(f,file_raw_sp4a,&ssp4a[0]);
    Write_Raw_Xmol(f,file_raw_sp4b,&ssp4b[0]);
    Write_Raw_Xmol(f,file_raw_6A,&s6A[0]);
    Write_Raw_Xmol(f,file_raw_6Z,&s6Z[0]);
    Write_Raw_Xmol(f,file_raw_7K,&s7K[0]);
    Write_Raw_Xmol(f,file_raw_sp5,&ssp5[0]);
    Write_Raw_Xmol(f,file_raw_sp5a,&ssp5a[0]);
    Write_Raw_Xmol(f,file_raw_sp5b,&ssp5b[0]);
    Write_Raw_Xmol(f,file_raw_7A,&s7A[0]);
    Write_Raw_Xmol(f,file_raw_8A,&s8A[0]);
    Write_Raw_Xmol(f,file_raw_8B,&s8B[0]);
    Write_Raw_Xmol(f,file_raw_8K,&s8K[0]);
    Write_Raw_Xmol(f,file_raw_9A,&s9A[0]);
    Write_Raw_Xmol(f,file_raw_9B,&s9B[0]);
    Write_Raw_Xmol(f,file_raw_9K,&s9K[0]);
    Write_Raw_Xmol(f,file_raw_10A,&s10A[0]);
    Write_Raw_Xmol(f,file_raw_10B,&s10B[0]);
    Write_Raw_Xmol(f,file_raw_10K,&s10K[0]);
    Write_Raw_Xmol(f,file_raw_10W,&s10W[0]);
    Write_Raw_Xmol(f,file_raw_11A,&s11A[0]);
    Write_Raw_Xmol(f,file_raw_11B,&s11B[0]);
    Write_Raw_Xmol(f,file_raw_11C,&s11C[0]);
    Write_Raw_Xmol(f,file_raw_11E,&s11E[0]);
    Write_Raw_Xmol(f,file_raw_11F,&s11F[0]);
    Write_Raw_Xmol(f,file_raw_11W,&s11W[0]);
    Write_Raw_Xmol(f,file_raw_12A,&s12A[0]);
    Write_Raw_Xmol(f,file_raw_12B,&s12B[0]);
    Write_Raw_Xmol(f,file_raw_12D,&s12D[0]);
    Write_Raw_Xmol(f,file_raw_12E,&s12E[0]);
    Write_Raw_Xmol(f,file_raw_12K,&s12K[0]);
    Write_Raw_Xmol(f,file_raw_13A,&s13A[0]);
    Write_Raw_Xmol(f,file_raw_13B,&s13B[0]);
    Write_Raw_Xmol(f,file_raw_13K,&s13K[0]);
    Write_Raw_Xmol(f,file_raw_FCC,&sFCC[0]);
    Write_Raw_Xmol(f,file_raw_HCP,&sHCP[0]);
    Write_Raw_Xmol(f,file_raw_BCC_9,&sBCC_9[0]);
    Write_Raw_Xmol(f,file_raw_BCC_15,&sBCC_15[0]);
}

void Write_Raw_Close() {
    fclose(file_raw_sp3);
    fclose(file_raw_sp3a);
    fclose(file_raw_sp3b);
    fclose(file_raw_5A);
    fclose(file_raw_sp4);
    fclose(file_raw_sp4a);
    fclose(file_raw_sp4b);
    fclose(file_raw_6A);
    fclose(file_raw_6Z);
    fclose(file_raw_7K);
    fclose(file_raw_sp5);
    fclose(file_raw_sp5a);
    fclose(file_raw_sp5b);
    fclose(file_raw_7A);
    fclose(file_raw_8A);
    fclose(file_raw_8B);
    fclose(file_raw_8K);
    fclose(file_raw_9A);
    fclose(file_raw_9B);
    fclose(file_raw_9K);
    fclose(file_raw_10A);
    fclose(file_raw_10B);
    fclose(file_raw_10K);
    fclose(file_raw_10W);
    fclose(file_raw_11A);
    fclose(file_raw_11B);
    fclose(file_raw_11C);
    fclose(file_raw_11E);
    fclose(file_raw_11F);
    fclose(file_raw_11W);
    fclose(file_raw_12A);
    fclose(file_raw_12B);
    fclose(file_raw_12D);
    fclose(file_raw_12E);
    fclose(file_raw_12K);
    fclose(file_raw_13A);
    fclose(file_raw_13B);
    fclose(file_raw_13K);
    fclose(file_raw_FCC);
    fclose(file_raw_HCP);
    fclose(file_raw_BCC_9);
    fclose(file_raw_BCC_15);
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

void Write_Cluster(int f, FILE *writeout, int *n, int **hc, int clusSize) {
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
    int i,j;
    int clusSize=5;

    fprintf(writeout,"%d\n",nsp5[f]);
    for (i=0;i<nsp5a[f];i++) {
        fprintf(writeout,"%d",sp5a[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5a[i][j]);
        fprintf(writeout,"	%d\n",sp5a[i][clusSize-1]);
    }
    for (i=0;i<nsp5b[f];i++) {
        fprintf(writeout,"%d",sp5b[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5b[i][j]);
        fprintf(writeout,"	%d\n",sp5b[i][clusSize-1]);
    }
    for (i=0;i<nsp5c[f];i++) {
        fprintf(writeout,"%d",sp5c[i][0]);
        for (j=1;j<clusSize-1;j++) fprintf(writeout,"	%d",sp5c[i][j]);
        fprintf(writeout,"	%d\n",sp5c[i][clusSize-1]);
    }
}