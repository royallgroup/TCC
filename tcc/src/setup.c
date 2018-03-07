#include "setup.h"
#include "tools.h"
#include "iniparser.h"
#include "globals.h"


void Setup_Output_Files() {
    char output_file[200];
    FILE *file_pointer;
    int cluster_number;

    if(do11AcenXyz == 1) {
        make_directory("centers_output");
        sprintf(output_file, "centers_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.%s_cen.xyz",
                fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[21]);
        file_pointer = open_file(output_file, "w");
        fclose(file_pointer);
    }

    if(do13AcenXyz == 1) {
        make_directory("centers_output");
        sprintf(output_file, "centers_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.%s_cen.xyz",
                fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[32]);
        file_pointer = open_file(output_file, "w");
        fclose(file_pointer);
    }

    if(doWriteBonds == 1) {
        sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
        file_pointer=fopen(output_file, "w");
        fclose(file_pointer);
    }

    if(doWritePopPerFrame == 1) {
        sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
        file_pointer=fopen(output_file, "w");
        fclose(file_pointer);
    }

    if(doWriteRaw == 1) {
        make_directory("raw_output");
        for(cluster_number=0; cluster_number < num_cluster_types; cluster_number++) {
            if (*do_cluster_list[cluster_number] == 1) {
                sprintf(output_file, "raw_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_%s",
                        fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
                file_pointer = open_file(output_file, "w");
                fclose(file_pointer);
            }
        }
    }

    if(doWriteClus == 1) {
        make_directory("cluster_output");
        for(cluster_number=0; cluster_number < num_cluster_types; cluster_number++) {
            if (*do_cluster_list[cluster_number] == 1) {
                sprintf(output_file, "cluster_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.clusts_%s",
                        fXmolName, rcutAA, rcutAB, rcutBB, Vor, fc, PBCs, cluster_names[cluster_number]);
                file_pointer = open_file(output_file, "w");
                fclose(file_pointer);
            }
        }
    }
}

void Initialise_Global_Variables() { // Initialize lots of important variables for static TCC algorithm
    int i, f, j;
    char errMsg[1000];

    dosp3=dosp3a=dosp3b=dosp3c=1;
    dosp4=dosp4a=dosp4b=dosp4c=1;
    dosp5=dosp5a=dosp5b=dosp5c=1;
    do6Z=do7K=do8A=do8B=do8K=do9A=do9B=do9K=do10A=do10B=do10K=do10W=1;
    do11A=do11B=do11C=do11E=do11F=do11W=do12A=do12B=do12D=do12E=do12K=1;
    do13A=do13B=do13K=doFCC=doHCP=doBCC9=1;
    doBCC15=0;
    
    msp3a=msp3b=msp3c=initNoStatic; // max size of **sp** arrays in dimension i
    msp4a=msp4b=msp4c=initNoStatic; // max size of **sp** arrays in dimension i
    msp5a=msp5b=msp5c=initNoStatic; // max size of **sp** arrays in dimension i
    m6Z=m7K=initNoStatic;   // max size of m** arrays in dimension i
    m8A=m8B=m8K=initNoStatic;   // max size of m** arrays in dimension i
    m9A=m9B=m9K=initNoStatic;   // max size of m** arrays in dimension i
    m10A=m10B=m10K=m10W=initNoStatic;   // max size of m** arrays in dimension i
    m11A=m11B=m11C=m11E=m11F=m11W=initNoStatic; // max size of m** arrays in dimension i
    m12A=m12B=m12D=m12E=m12K=initNoStatic;  // max size of m** arrays in dimension i
    m13A=m13B=m13K=initNoStatic;    // max size of m** arrays in dimension i
    mFCC=mHCP=mBCC_9=mBCC_15=initNoStatic;  // max size of **sp** arrays in dimension i

    mmem_sp3b=mmem_sp3c=mmem_sp4b=mmem_sp4c=mmem_sp5b=mmem_sp5c=initNoClustPerPart;
    
    mean_pop_per_frame = malloc(num_cluster_types*sizeof(double));
    for(i=0; i<num_cluster_types; i++) {
        mean_pop_per_frame[i] = 0.0;
    }

    tiltxy = tiltxz = tiltyz = 0;
    x = malloc(max_particle_number*sizeof(double));   if (x==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): x[] malloc out of memory\n");    Error_no_free(errMsg); }    // positions of particles in a configuration
    y = malloc(max_particle_number*sizeof(double));   if (y==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): y[] malloc out of memory\n");    Error_no_free(errMsg); }
    z = malloc(max_particle_number*sizeof(double));   if (z==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): z[] malloc out of memory\n");    Error_no_free(errMsg); }
    particle_type=malloc(max_particle_number*sizeof(int)); if (particle_type==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): particle_type[] malloc out of memory\n");   Error_no_free(errMsg); }    // type of species

    num_bonds = malloc(max_particle_number*sizeof(int));    if (num_bonds==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): num_bonds[] malloc out of memory\n");    Error_no_free(errMsg); }    // number of "bonded" neighbours of a particle

    bNums = malloc(max_particle_number*sizeof(int *));    if (bNums==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): bNums[] malloc out of memory\n");    Error_no_free(errMsg); }    // list of bonded particles to each particle
    for (j=0; j<max_particle_number; ++j) { bNums[j] = malloc(nB*sizeof(int));    if (bNums[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): bNums[][] malloc out of memory\n");   Error_no_free(errMsg); } }

    squared_bondlengths = malloc(max_particle_number*sizeof(double *));   if (squared_bondlengths==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): squared_bondlengths[] malloc out of memory\n");    Error_no_free(errMsg); }    // array of bond lengths
    for (j=0; j<max_particle_number; ++j) { squared_bondlengths[j] = malloc(nB*sizeof(double));   if (squared_bondlengths[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): squared_bondlengths[][] malloc out of memory\n");   Error_no_free(errMsg); } }

    // arrays for the number of clusters of each type bonded to each particle
    mem_sp3b = malloc(max_particle_number*sizeof(int *)); if (mem_sp3b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp3b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp3b[j] = malloc(mmem_sp3b*sizeof(int));  if (mem_sp3b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp3b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp3b = malloc(max_particle_number*sizeof(int));  if (nmem_sp3b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp3b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp3c = malloc(max_particle_number*sizeof(int *)); if (mem_sp3c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp3c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp3c[j] = malloc(mmem_sp3c*sizeof(int));  if (mem_sp3c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp3c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp3c = malloc(max_particle_number*sizeof(int));  if (nmem_sp3c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp3c[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp4b = malloc(max_particle_number*sizeof(int *)); if (mem_sp4b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp4b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp4b[j] = malloc(mmem_sp4b*sizeof(int));  if (mem_sp4b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp4b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp4b = malloc(max_particle_number*sizeof(int));  if (nmem_sp4b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp4b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp4c = malloc(max_particle_number*sizeof(int *)); if (mem_sp4c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp4c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp4c[j] = malloc(mmem_sp4c*sizeof(int));  if (mem_sp4c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp4c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp4c = malloc(max_particle_number*sizeof(int));  if (nmem_sp4c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp4c[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp5b = malloc(max_particle_number*sizeof(int *)); if (mem_sp5b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp5b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp5b[j] = malloc(mmem_sp5b*sizeof(int));  if (mem_sp5b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp5b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp5b = malloc(max_particle_number*sizeof(int));  if (nmem_sp5b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp5b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp5c = malloc(max_particle_number*sizeof(int *)); if (mem_sp5c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp5c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<max_particle_number; ++j) { mem_sp5c[j] = malloc(mmem_sp5c*sizeof(int));  if (mem_sp5c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): mem_sp5c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp5c = malloc(max_particle_number*sizeof(int));  if (nmem_sp5c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): nmem_sp5c[] malloc out of memory\n");    Error_no_free(errMsg); }

    // arrays for the clusters found in each frame
    hcsp3a = malloc(msp3a*sizeof(int *)); if (hcsp3a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3a; ++j) { hcsp3a[j] = malloc(3*sizeof(int));  if (hcsp3a[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp3b = malloc(msp3b*sizeof(int *)); if (hcsp3b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3b; ++j) { hcsp3b[j] = malloc(4*sizeof(int));  if (hcsp3b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp3c = malloc(msp3c*sizeof(int *)); if (hcsp3c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3c; ++j) { hcsp3c[j] = malloc(5*sizeof(int));  if (hcsp3c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp3c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp4a = malloc(msp4a*sizeof(int *)); if (hcsp4a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4a; ++j) { hcsp4a[j] = malloc(4*sizeof(int));  if (hcsp4a[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp4b = malloc(msp4b*sizeof(int *)); if (hcsp4b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4b; ++j) { hcsp4b[j] = malloc(5*sizeof(int));  if (hcsp4b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp4c = malloc(msp4c*sizeof(int *)); if (hcsp4c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4c; ++j) { hcsp4c[j] = malloc(6*sizeof(int));  if (hcsp4c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp4c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp5a = malloc(msp5a*sizeof(int *)); if (hcsp5a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5a; ++j) { hcsp5a[j] = malloc(5*sizeof(int));  if (hcsp5a[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp5b = malloc(msp5b*sizeof(int *)); if (hcsp5b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5b; ++j) { hcsp5b[j] = malloc(6*sizeof(int));  if (hcsp5b[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hcsp5c = malloc(msp5c*sizeof(int *)); if (hcsp5c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5c; ++j) { hcsp5c[j] = malloc(7*sizeof(int));  if (hcsp5c[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sp5c[][] malloc out of memory\n"); Error_no_free(errMsg); } }


    // arrays for the clusters found in each frame
    hc6Z = malloc(m6Z*sizeof(int *));   if (hc6Z==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc6Z[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m6Z; ++j) { hc6Z[j] = malloc(6*sizeof(int));    if (hc6Z[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc6Z[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc7K = malloc(m7K*sizeof(int *));   if (hc7K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc7K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m7K; ++j) { hc7K[j] = malloc(7*sizeof(int));    if (hc7K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc7K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc8A = malloc(m8A*sizeof(int *));   if (hc8A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8A[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8A; ++j) { hc8A[j] = malloc(8*sizeof(int));    if (hc8A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8A[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc8B = malloc(m8B*sizeof(int *));   if (hc8B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8B[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8B; ++j) { hc8B[j] = malloc(8*sizeof(int));    if (hc8B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8B[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc8K = malloc(m8K*sizeof(int *));   if (hc8K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8K; ++j) { hc8K[j] = malloc(8*sizeof(int));    if (hc8K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc8K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc9A = malloc(m9A*sizeof(int *));   if (hc9A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9A[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9A; ++j) { hc9A[j] = malloc(9*sizeof(int));    if (hc9A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9A[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc9B = malloc(m9B*sizeof(int *));   if (hc9B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9B[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9B; ++j) { hc9B[j] = malloc(9*sizeof(int));    if (hc9B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9B[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc9K = malloc(m9K*sizeof(int *));   if (hc9K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9K; ++j) { hc9K[j] = malloc(9*sizeof(int));    if (hc9K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc9K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc10A = malloc(m10A*sizeof(int *)); if (hc10A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10A; ++j) { hc10A[j] = malloc(10*sizeof(int)); if (hc10A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10B = malloc(m10B*sizeof(int *)); if (hc10B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10B; ++j) { hc10B[j] = malloc(10*sizeof(int)); if (hc10B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10K = malloc(m10K*sizeof(int *)); if (hc10K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10K; ++j) { hc10K[j] = malloc(10*sizeof(int)); if (hc10K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10W = malloc(m10W*sizeof(int *)); if (hc10W==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10W[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10W; ++j) { hc10W[j] = malloc(10*sizeof(int)); if (hc10W[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc10W[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hc11A = malloc(m11A*sizeof(int *)); if (hc11A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11A; ++j) { hc11A[j] = malloc(11*sizeof(int)); if (hc11A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11B = malloc(m11B*sizeof(int *)); if (hc11B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11B; ++j) { hc11B[j] = malloc(11*sizeof(int)); if (hc11B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11C = malloc(m11C*sizeof(int *)); if (hc11C==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11C[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11C; ++j) { hc11C[j] = malloc(11*sizeof(int)); if (hc11C[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11C[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11E = malloc(m11E*sizeof(int *)); if (hc11E==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11E[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11E; ++j) { hc11E[j] = malloc(11*sizeof(int)); if (hc11E[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11E[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11F = malloc(m11F*sizeof(int *)); if (hc11F==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11F[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11F; ++j) { hc11F[j] = malloc(11*sizeof(int)); if (hc11F[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11F[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11W = malloc(m11W*sizeof(int *)); if (hc11W==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11W[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11W; ++j) { hc11W[j] = malloc(11*sizeof(int)); if (hc11W[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc11W[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hc12A = malloc(m12A*sizeof(int *)); if (hc12A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12A; ++j) { hc12A[j] = malloc(12*sizeof(int)); if (hc12A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12B = malloc(m12B*sizeof(int *)); if (hc12B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12B; ++j) { hc12B[j] = malloc(12*sizeof(int)); if (hc12B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12D = malloc(m12D*sizeof(int *)); if (hc12D==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12D[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12D; ++j) { hc12D[j] = malloc(12*sizeof(int)); if (hc12D[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12D[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12E = malloc(m12E*sizeof(int *)); if (hc12E==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12E[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12E; ++j) { hc12E[j] = malloc(12*sizeof(int)); if (hc12E[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12E[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12K = malloc(m12K*sizeof(int *)); if (hc12K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12K; ++j) { hc12K[j] = malloc(12*sizeof(int)); if (hc12K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc12K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
        
    hc13A = malloc(m13A*sizeof(int *)); if (hc13A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13A; ++j) { hc13A[j] = malloc(13*sizeof(int)); if (hc13A[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc13B = malloc(m13B*sizeof(int *)); if (hc13B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13B; ++j) { hc13B[j] = malloc(13*sizeof(int)); if (hc13B[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc13K = malloc(m13K*sizeof(int *)); if (hc13K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13K; ++j) { hc13K[j] = malloc(13*sizeof(int)); if (hc13K[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hc13K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hcFCC = malloc(mFCC*sizeof(int *)); if (hcFCC==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcFCC[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mFCC; ++j) { hcFCC[j] = malloc(13*sizeof(int)); if (hcFCC[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcFCC[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcHCP = malloc(mHCP*sizeof(int *)); if (hcHCP==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcHCP[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mHCP; ++j) { hcHCP[j] = malloc(13*sizeof(int)); if (hcHCP[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcHCP[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcBCC_9 = malloc(mBCC_9*sizeof(int *)); if (hcBCC_9==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcBCC_9[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mBCC_9; ++j) { hcBCC_9[j] = malloc(9*sizeof(int));  if (hcBCC_9[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcBCC_9[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcBCC_15 = malloc(mBCC_15*sizeof(int *));   if (hcBCC_15==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcBCC_15[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<mBCC_15; ++j) { hcBCC_15[j] = malloc(15*sizeof(int));   if (hcBCC_15[j]==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): hcBCC_15[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    // character arrays listing what type each particle is when found in a cluster
    ssp3a=malloc(max_particle_number*sizeof(char)); if (ssp3a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp3a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp3b=malloc(max_particle_number*sizeof(char)); if (ssp3b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp3b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp3c=malloc(max_particle_number*sizeof(char)); if (ssp3c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp3c[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp4a=malloc(max_particle_number*sizeof(char)); if (ssp4a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp4a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4b=malloc(max_particle_number*sizeof(char)); if (ssp4b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp4b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4c=malloc(max_particle_number*sizeof(char)); if (ssp4c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp4c[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp5a=malloc(max_particle_number*sizeof(char)); if (ssp5a==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp5a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5b=malloc(max_particle_number*sizeof(char)); if (ssp5b==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp5b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5c=malloc(max_particle_number*sizeof(char)); if (ssp5c==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): ssp5c[] malloc out of memory\n"); Error_no_free(errMsg); }

    s6Z=malloc(max_particle_number*sizeof(char)); if (s6Z==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s6Z[] malloc out of memory\n"); Error_no_free(errMsg); }
    s7K=malloc(max_particle_number*sizeof(char)); if (s7K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s7K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s8A=malloc(max_particle_number*sizeof(char)); if (s8A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s8A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s8B=malloc(max_particle_number*sizeof(char)); if (s8B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s8B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s8K=malloc(max_particle_number*sizeof(char)); if (s8K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s8K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9A=malloc(max_particle_number*sizeof(char)); if (s9A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s9A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9B=malloc(max_particle_number*sizeof(char)); if (s9B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s9B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9K=malloc(max_particle_number*sizeof(char)); if (s9K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s9K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10A=malloc(max_particle_number*sizeof(char)); if (s10A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s10A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10B=malloc(max_particle_number*sizeof(char)); if (s10B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s10B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10K=malloc(max_particle_number*sizeof(char)); if (s10K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s10K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10W=malloc(max_particle_number*sizeof(char)); if (s10W==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s10W[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11A=malloc(max_particle_number*sizeof(char)); if (s11A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11B=malloc(max_particle_number*sizeof(char)); if (s11B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11C=malloc(max_particle_number*sizeof(char)); if (s11C==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11C[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11E=malloc(max_particle_number*sizeof(char)); if (s11E==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11E[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11F=malloc(max_particle_number*sizeof(char)); if (s11F==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11F[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11W=malloc(max_particle_number*sizeof(char)); if (s11W==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s11W[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12A=malloc(max_particle_number*sizeof(char)); if (s12A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s12A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12B=malloc(max_particle_number*sizeof(char)); if (s12B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s12B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12D=malloc(max_particle_number*sizeof(char)); if (s12D==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s12D[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12E=malloc(max_particle_number*sizeof(char)); if (s12E==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s12E[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12K=malloc(max_particle_number*sizeof(char)); if (s12K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s12K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13A=malloc(max_particle_number*sizeof(char)); if (s13A==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s13A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13B=malloc(max_particle_number*sizeof(char)); if (s13B==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s13B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13K=malloc(max_particle_number*sizeof(char)); if (s13K==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): s13K[] malloc out of memory\n"); Error_no_free(errMsg); }
    sFCC=malloc(max_particle_number*sizeof(char)); if (sFCC==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sFCC[] malloc out of memory\n"); Error_no_free(errMsg); }
    sHCP=malloc(max_particle_number*sizeof(char)); if (sHCP==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sHCP[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_9=malloc(max_particle_number*sizeof(char)); if (sBCC_9==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sBCC_9[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_15=malloc(max_particle_number*sizeof(char)); if (sBCC_15==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): sBCC_15[] malloc out of memory\n"); Error_no_free(errMsg); }

    // particle fraction of particles in each cluster in each frame
    pop_per_frame = malloc(num_cluster_types*sizeof(double *));

    for(i=0; i<num_cluster_types; i++) {
        pop_per_frame[i] = malloc(FRAMES*sizeof(double));
        if (pop_per_frame[i]==NULL) {
            sprintf(errMsg,"Initialise_Global_Variables(): pop_per_frame malloc out of memory\n");
            Error_no_free(errMsg);
        }
    }

    // Zero the pop per frame arrays
    for(i=0; i<num_cluster_types; i++) {
        for (f=0; f<FRAMES; f++) {
            pop_per_frame[i][f] = 0;
        }
    }
}

void Reset_Frame_Variables() { // Reset static variables in each frame
    int i, cluster_type;

    for(cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        *num_cluster_list[cluster_type] = 0;
    }

    memset(num_bonds, 0, current_frame_particle_number* sizeof(int));

    memset(nmem_sp3b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp3c, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp4b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp4c, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp5b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp5c, 0, sizeof(int)*max_particle_number);

    for(cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        for (i = 0; i < max_particle_number; ++i) {
            memset(*raw_list[cluster_type], 'C', max_particle_number*sizeof(char));
        }
    }
}

void Free_All_Variables()  {  // Free bond detection variables
    int i;

    free(mean_pop_per_frame);

    free(particle_type);
    free(fXmolName);
    free(fBoxSizeName);
    free(x); free(y); free(z);

    for (i=0; i<max_particle_number; ++i) {
        free(bNums[i]); 
        free(squared_bondlengths[i]);
    }
    free(bNums); free(squared_bondlengths); free(num_bonds);
    
    for (i=0; i<num_cluster_types; ++i) {
        free(pop_per_frame[i]);
    }
    free(pop_per_frame);

    for (i=0; i<msp3a; ++i) free(hcsp3a[i]);
    for (i=0; i<msp3b; ++i) free(hcsp3b[i]);
    for (i=0; i<msp3c; ++i) free(hcsp3c[i]);
    for (i=0; i<msp4a; ++i) free(hcsp4a[i]);
    for (i=0; i<msp4b; ++i) free(hcsp4b[i]);
    for (i=0; i<msp4c; ++i) free(hcsp4c[i]);
    for (i=0; i<msp5a; ++i) free(hcsp5a[i]);
    for (i=0; i<msp5b; ++i) free(hcsp5b[i]);
    for (i=0; i<msp5c; ++i) free(hcsp5c[i]);
    for (i=0; i<m6Z; ++i) free(hc6Z[i]);
    for (i=0; i<m7K; ++i) free(hc7K[i]);
    for (i=0; i<m8A; ++i) free(hc8A[i]);
    for (i=0; i<m8B; ++i) free(hc8B[i]);
    for (i=0; i<m8K; ++i) free(hc8K[i]);
    for (i=0; i<m9A; ++i) free(hc9A[i]);
    for (i=0; i<m9B; ++i) free(hc9B[i]);
    for (i=0; i<m9K; ++i) free(hc9K[i]);
    for (i=0; i<m10A; ++i) free(hc10A[i]);
    for (i=0; i<m10B; ++i) free(hc10B[i]);
    for (i=0; i<m10K; ++i) free(hc10K[i]);
    for (i=0; i<m10W; ++i) free(hc10W[i]);
    for (i=0; i<m11A; ++i) free(hc11A[i]);
    for (i=0; i<m11B; ++i) free(hc11B[i]);
    for (i=0; i<m11C; ++i) free(hc11C[i]);
    for (i=0; i<m11E; ++i) free(hc11E[i]);
    for (i=0; i<m11F; ++i) free(hc11F[i]);
    for (i=0; i<m11W; ++i) free(hc11W[i]);
    for (i=0; i<m12A; ++i) free(hc12A[i]);
    for (i=0; i<m12B; ++i) free(hc12B[i]);
    for (i=0; i<m12D; ++i) free(hc12D[i]);
    for (i=0; i<m12E; ++i) free(hc12E[i]);
    for (i=0; i<m12K; ++i) free(hc12K[i]);
    for (i=0; i<m13A; ++i) free(hc13A[i]);
    for (i=0; i<m13B; ++i) free(hc13B[i]);
    for (i=0; i<m13K; ++i) free(hc13K[i]);
    for (i=0; i<mFCC; ++i) free(hcFCC[i]);
    for (i=0; i<mHCP; ++i) free(hcHCP[i]);
    for (i=0; i<mBCC_9; ++i) free(hcBCC_9[i]);
    for (i=0; i<mBCC_15; ++i) free(hcBCC_15[i]);
    
    for (i=0; i<max_particle_number; ++i) free(mem_sp3b[i]);
    for (i=0; i<max_particle_number; ++i) free(mem_sp3c[i]);
    for (i=0; i<max_particle_number; ++i) free(mem_sp4b[i]);
    for (i=0; i<max_particle_number; ++i) free(mem_sp4c[i]);
    for (i=0; i<max_particle_number; ++i) free(mem_sp5b[i]);
    for (i=0; i<max_particle_number; ++i) free(mem_sp5c[i]);
    
    free(mem_sp3b);
    free(mem_sp3c);
    free(mem_sp4b);
    free(mem_sp4c);
    free(mem_sp5b);
    free(mem_sp5c);

    free(nmem_sp3b);
    free(nmem_sp3c);
    free(nmem_sp4b);
    free(nmem_sp4c);
    free(nmem_sp5b);
    free(nmem_sp5c);

    free(hcsp3a); free(hcsp3b); free(hcsp3c);
    free(hcsp4a); free(hcsp4b); free(hcsp4c);
    free(hcsp5a); free(hcsp5b); free(hcsp5c);
    free(hc6Z); free(hc7K);
    free(hc8A); free(hc8B); free(hc8K);
    free(hc9A); free(hc9B); free(hc9K);
    free(hc10A); free(hc10B); free(hc10K); free(hc10W);
    free(hc11A); free(hc11B); free(hc11C); free(hc11E); free(hc11F); free(hc11W);
    free(hc12A); free(hc12B); free(hc12D); free(hc12E); free(hc12K);
    free(hc13A); free(hc13B); free(hc13K);
    free(hcFCC); free(hcHCP); free(hcBCC_9); free(hcBCC_15);

    free(ssp3a); free(ssp3b); free(ssp3c);
    free(ssp4a); free(ssp4b); free(ssp4c);
    free(s6Z); free(s7K);
    free(ssp5a); free(ssp5b); free(ssp5c);
    free(s8A); free(s8B); free(s8K);
    free(s9A); free(s9B); free(s9K);
    free(s10A); free(s10B); free(s10K); free(s10W);
    free(s11A); free(s11B); free(s11C); free(s11E); free(s11F); free(s11W);
    free(s12A); free(s12B); free(s12D); free(s12E); free(s12K);
    free(s13A); free(s13B); free(s13K);
    free(sFCC);free(sHCP); free(sBCC_9); free(sBCC_15);
}
