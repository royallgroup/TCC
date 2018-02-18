#include "stats.h"
#include "globals.h"
#include "tools.h"

void Stats_Init() {
    int i;
    char errMsg[1000];

    a5=malloc(N*sizeof(int));	if (a5==NULL) { sprintf(errMsg,"Stats_Init(): a5[] malloc out of memory\n");	Error(errMsg); }
    a6=malloc(N*sizeof(int));	if (a6==NULL) { sprintf(errMsg,"Stats_Init(): a6[] malloc out of memory\n");	Error(errMsg); }
    a7=malloc(N*sizeof(int));	if (a7==NULL) { sprintf(errMsg,"Stats_Init(): a7[] malloc out of memory\n");	Error(errMsg); }
    a8=malloc(N*sizeof(int));	if (a8==NULL) { sprintf(errMsg,"Stats_Init(): a8[] malloc out of memory\n");	Error(errMsg); }
    a9=malloc(N*sizeof(int));	if (a9==NULL) { sprintf(errMsg,"Stats_Init(): a9[] malloc out of memory\n");	Error(errMsg); }
    a10=malloc(N*sizeof(int));	if (a10==NULL) { sprintf(errMsg,"Stats_Init(): a10[] malloc out of memory\n");	Error(errMsg); }
    a11=malloc(N*sizeof(int));	if (a11==NULL) { sprintf(errMsg,"Stats_Init(): a11[] malloc out of memory\n");	Error(errMsg); }
    a12=malloc(N*sizeof(int));	if (a12==NULL) { sprintf(errMsg,"Stats_Init(): a12[] malloc out of memory\n");	Error(errMsg); }
    a13=malloc(N*sizeof(int));	if (a13==NULL) { sprintf(errMsg,"Stats_Init(): a13[] malloc out of memory\n");	Error(errMsg); }
    a15=malloc(N*sizeof(int));	if (a15==NULL) { sprintf(errMsg,"Stats_Init(): a15[] malloc out of memory\n");	Error(errMsg); }

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }

    net_clusters = malloc(num_cluster_types*sizeof(int));
    gross_clusters = malloc(num_cluster_types*sizeof(int));
    total_clusters = malloc(num_cluster_types*sizeof(int));

    for(i=0; i<num_cluster_types; i++) {
        net_clusters[i] = 0;
        gross_clusters[i] = 0;
        total_clusters[i] = 0;
    }
}

void Stats_Reset() {
    int i;

    for(i=0;i<N;i++) {
        a5[i]=a6[i]=a7[i]=a8[i]=a9[i]=a10[i]=a11[i]=a12[i]=a13[i]=a15[i]=0;
    }
}

void Stats_FreeMem() {
    free(a5);
    free(a6);
    free(a7);
    free(a8);
    free(a9);
    free(a10);
    free(a11);
    free(a12);
    free(a13);
    free(a15);
    free(net_clusters);
    free(gross_clusters);
    free(total_clusters);
}

void Stats_SetA() { // Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles
    int flg1, flg2;
    int i;

    for (i=0;i<N;i++) {
        flg1 = sFCC[i] != 'C' || sHCP[i] != 'C' || sBCC_15[i] != 'C' || s13A[i] != 'C' || s13B[i] != 'C' || s13K[i] != 'C';
        if(flg1==1) a13[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s12A[i] != 'C' || s12B[i] != 'C' || s12D[i] != 'C' || s12E[i] != 'C' || s12K[i] != 'C';
        flg2 = a13[i];
        if(flg1==1 || flg2==1) a12[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s11A[i] != 'C' || s11B[i] != 'C' || s11C[i] != 'C' || s11E[i] != 'C' || s11F[i] != 'C' || s11W[i] != 'C';
        flg2 = a12[i];
        if(flg1==1 || flg2==1) a11[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s10A[i] != 'C' || s10B[i] != 'C' || s10K[i] != 'C' || s10W[i] != 'C';
        flg2 = a11[i];
        if(flg1==1 || flg2==1) a10[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s9A[i] != 'C' || s9B[i] != 'C' || s9K[i] != 'C' || sBCC_9[i] != 'C';
        flg2 = a10[i];
        if(flg1==1 || flg2==1) a9[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = s8A[i] != 'C' || s8B[i] != 'C' || s8K[i] != 'C';
        flg2 = a9[i];
        if(flg1==1 || flg2==1) a8[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = ssp5c[i] != 'C' || s7K[i] != 'C';
        flg2 = a8[i];
        if(flg1==1 || flg2==1) a7[i] = 1;
    }
    for (i=0;i<N;i++) {
        flg1 = ssp4c[i] != 'C' || s6Z[i] != 'C';
        flg2 = a7[i];
        if(flg1==1 || flg2==1) a6[i] = 1;
    }
}

void Stats_Analyse() {
    int i;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') ++gross_clusters[0];
        if(ssp3a[i] != 'C') ++gross_clusters[1];
        if(ssp3b[i] != 'C') ++gross_clusters[2];
        if(ssp3c[i] != 'C') ++gross_clusters[3];
        if(ssp4[i] != 'C') ++gross_clusters[4];
        if(ssp4a[i] != 'C') ++gross_clusters[5];
        if(ssp4b[i] != 'C') ++gross_clusters[6];
        if(ssp4c[i] != 'C') ++gross_clusters[7];
        if(s6Z[i] != 'C') ++gross_clusters[8];
        if(s7K[i] != 'C') ++gross_clusters[9];
        if(ssp5[i] != 'C') ++gross_clusters[10];
        if(ssp5a[i] != 'C') ++gross_clusters[11];
        if(ssp5b[i] != 'C') ++gross_clusters[12];
        if(ssp5c[i] != 'C') ++gross_clusters[13];
        if(s8A[i] != 'C') ++gross_clusters[14];
        if(s8B[i] != 'C') ++gross_clusters[15];
        if(s8K[i] != 'C') ++gross_clusters[16];
        if(s9A[i] != 'C') ++gross_clusters[17];
        if(s9B[i] != 'C') ++gross_clusters[18];
        if(s9K[i] != 'C') ++gross_clusters[19];
        if(s10A[i] != 'C') ++gross_clusters[20];
        if(s10B[i] != 'C') ++gross_clusters[21];
        if(s10K[i] != 'C') ++gross_clusters[22];
        if(s10W[i] != 'C') ++gross_clusters[23];
        if(s11A[i] != 'C') ++gross_clusters[24];
        if(s11B[i] != 'C') ++gross_clusters[25];
        if(s11C[i] != 'C') ++gross_clusters[26];
        if(s11E[i] != 'C') ++gross_clusters[27];
        if(s11F[i] != 'C') ++gross_clusters[28];
        if(s11W[i] != 'C') ++gross_clusters[29];
        if(s12A[i] != 'C') ++gross_clusters[30];
        if(s12B[i] != 'C') ++gross_clusters[31];
        if(s12D[i] != 'C') ++gross_clusters[32];
        if(s12E[i] != 'C') ++gross_clusters[33];
        if(s12K[i] != 'C') ++gross_clusters[34];
        if(s13A[i] != 'C') ++gross_clusters[35];
        if(s13B[i] != 'C') ++gross_clusters[36];
        if(s13K[i] != 'C') ++gross_clusters[37];
        if(sFCC[i] != 'C') ++gross_clusters[38];
        if(sHCP[i] != 'C') ++gross_clusters[39];
        if(sBCC_9[i] != 'C') ++gross_clusters[40];
        if(sBCC_15[i] != 'C') ++gross_clusters[41];
    }
    Stats_SetA();
    for(i=0; i<N; ++i){
        if(ssp3c[i] != 'C' && !a6[i]) ++net_clusters[3];
        if(ssp4c[i] != 'C' && !a7[i]) ++net_clusters[7];
        if(s6Z[i] != 'C' && !a7[i]) ++net_clusters[8];
        if(s7K[i] != 'C' && !a7[i]) ++net_clusters[9];
        if(ssp5c[i] != 'C' && !a8[i]) ++net_clusters[13];
        if(s8A[i] != 'C' && !a9[i]) ++net_clusters[14];
        if(s8B[i] != 'C' && !a9[i]) ++net_clusters[15];
        if(s8K[i] != 'C' && !a9[i]) ++net_clusters[16];
        if(s9A[i] != 'C' && !a10[i]) ++net_clusters[17];
        if(s9B[i] != 'C' && !a10[i]) ++net_clusters[18];
        if(s9K[i] != 'C' && !a10[i]) ++net_clusters[19];
        if(s10A[i] != 'C' && !a11[i]) ++net_clusters[20];
        if(s10B[i] != 'C' && !a11[i]) ++net_clusters[21];
        if(s10K[i] != 'C' && !a11[i]) ++net_clusters[22];
        if(s10W[i] != 'C' && !a11[i]) ++net_clusters[23];
        if(s11A[i] != 'C' && !a12[i]) ++net_clusters[24];
        if(s11B[i] != 'C' && !a12[i]) ++net_clusters[25];
        if(s11C[i] != 'C' && !a12[i]) ++net_clusters[26];
        if(s11E[i] != 'C' && !a12[i]) ++net_clusters[27];
        if(s11F[i] != 'C' && !a12[i]) ++net_clusters[28];
        if(s11W[i] != 'C' && !a12[i]) ++net_clusters[29];
        if(s12A[i] != 'C' && !a13[i]) ++net_clusters[30];
        if(s12B[i] != 'C' && !a13[i]) ++net_clusters[31];
        if(s12D[i] != 'C' && !a13[i]) ++net_clusters[32];
        if(s12E[i] != 'C' && !a13[i]) ++net_clusters[33];
        if(s12K[i] != 'C' && !a13[i]) ++net_clusters[34];
        if(s13A[i] != 'C') ++net_clusters[35];
        if(s13B[i] != 'C') ++net_clusters[36];
        if(s13K[i] != 'C') ++net_clusters[37];
        if(sFCC[i] != 'C') ++net_clusters[38];
        if(sHCP[i] != 'C') ++net_clusters[39];
        if(sBCC_9[i] != 'C') ++net_clusters[40];
        if(sBCC_15[i] != 'C') ++net_clusters[41];
    }
}

void Accuumlate_Stats() {
    int i;

    for(i=0; i<num_cluster_types; i++) {
        total_clusters[i] += *num_cluster_list[i];
    }
}

void Stats_Report(char *filename) {
    char errMsg[1000];
    FILE *writeout;
    int i;

    printf("Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    for(i=0; i<num_cluster_types; i++) {
        printf("%s	%d	%d	%d	%.5lg\n",
               cluster_names[i], total_clusters[i], gross_clusters[i], net_clusters[i], mean_pop_per_frame[i]);
    }
    printf("maxnB	%d\n",maxnb);
    printf("correctedBonds %d per frame %.5lg per part per frame %.5lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));


    writeout=fopen(filename,"w");
    if (writeout==NULL)  {
        sprintf(errMsg,"Stats_Report(): Error opening file %s",filename);	// Always test file open
        Error(errMsg);
    }
    fprintf(writeout,"%s\n",filename);
    fprintf(writeout,"Clust	Number	Gross	Net	Mean Pop Per Frame\n");
    for(i=0; i<num_cluster_types; i++) {
        fprintf(writeout, "%s	%d	%d	%d	%.5lg\n", cluster_names[i],
                total_clusters[i], gross_clusters[i], net_clusters[i], mean_pop_per_frame[i]);
    }
    fprintf(writeout,"maxnB	%d\n",maxnb);
    fprintf(writeout,"correctedBonds	%d	per frame	%.15lg	per part per frame	%.15lg\n",correctedBonds,(double)correctedBonds/FRAMES,(double)correctedBonds/(FRAMES*N));

    fclose(writeout);
}

void Pop_Per_Frame(int f) {
    int i;

    for(i=0; i<N; ++i){
        if(ssp3[i] != 'C') pop_per_frame[0][f]+=1.0;
        if(ssp3a[i] != 'C') pop_per_frame[1][f]+=1.0;
        if(ssp3b[i] != 'C') pop_per_frame[2][f]+=1.0;
        if(ssp3c[i] != 'C') pop_per_frame[3][f]+=1.0;
        if(ssp4[i] != 'C') pop_per_frame[4][f]+=1.0;
        if(ssp4a[i] != 'C') pop_per_frame[5][f]+=1.0;
        if(ssp4b[i] != 'C') pop_per_frame[6][f]+=1.0;
        if(ssp4c[i] != 'C') pop_per_frame[7][f]+=1.0;
        if(ssp5[i] != 'C') pop_per_frame[8][f]+=1.0;
        if(ssp5a[i] != 'C') pop_per_frame[9][f]+=1.0;
        if(ssp5b[i] != 'C') pop_per_frame[10][f]+=1.0;
        if(ssp5c[i] != 'C') pop_per_frame[11][f]+=1.0;
        if(s6Z[i] != 'C') pop_per_frame[12][f]+=1.0;
        if(s7K[i] != 'C') pop_per_frame[13][f]+=1.0;
        if(s8A[i] != 'C') pop_per_frame[14][f]+=1.0;
        if(s8B[i] != 'C') pop_per_frame[15][f]+=1.0;
        if(s8K[i] != 'C') pop_per_frame[16][f]+=1.0;
        if(s9A[i] != 'C') pop_per_frame[17][f]+=1.0;
        if(s9B[i] != 'C') pop_per_frame[18][f]+=1.0;
        if(s9K[i] != 'C') pop_per_frame[19][f]+=1.0;
        if(s10A[i] != 'C') pop_per_frame[20][f]+=1.0;
        if(s10B[i] != 'C') pop_per_frame[21][f]+=1.0;
        if(s10K[i] != 'C') pop_per_frame[22][f]+=1.0;
        if(s10W[i] != 'C') pop_per_frame[23][f]+=1.0;
        if(s11A[i] != 'C') pop_per_frame[24][f]+=1.0;
        if(s11B[i] != 'C') pop_per_frame[25][f]+=1.0;
        if(s11C[i] != 'C') pop_per_frame[26][f]+=1.0;
        if(s11E[i] != 'C') pop_per_frame[27][f]+=1.0;
        if(s11F[i] != 'C') pop_per_frame[28][f]+=1.0;
        if(s11W[i] != 'C') pop_per_frame[29][f]+=1.0;
        if(s12A[i] != 'C') pop_per_frame[30][f]+=1.0;
        if(s12B[i] != 'C') pop_per_frame[31][f]+=1.0;
        if(s12D[i] != 'C') pop_per_frame[32][f]+=1.0;
        if(s12E[i] != 'C') pop_per_frame[33][f]+=1.0;
        if(s12K[i] != 'C') pop_per_frame[34][f]+=1.0;
        if(s13A[i] != 'C') pop_per_frame[35][f]+=1.0;
        if(s13B[i] != 'C') pop_per_frame[36][f]+=1.0;
        if(s13K[i] != 'C') pop_per_frame[37][f]+=1.0;
        if(sFCC[i] != 'C') pop_per_frame[38][f]+=1.0;
        if(sHCP[i] != 'C') pop_per_frame[39][f]+=1.0;
        if(sBCC_9[i] != 'C') pop_per_frame[40][f]+=1.0;
        if(sBCC_15[i] != 'C') pop_per_frame[41][f]+=1.0;
    }

    // Add pop per frame to running total over all frames
    for(i=0; i<num_cluster_types; i++) {
        pop_per_frame[i][f]/=N;
        mean_pop_per_frame[i] += pop_per_frame[i][f];
    }
}

void Normalise_Populations() {
    int i;
    for(i=0; i<num_cluster_types; i++) {
        mean_pop_per_frame[i]/=FRAMES;
    }
}