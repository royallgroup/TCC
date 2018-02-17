#include "setup.h"
#include "tools.h"
#include "iniparser.h"
#include "globals.h"
#include "math.h"

void Setup_ReadIniFile(char *filename) {
    
    char errMsg[1000];
    double RHO;
    dictionary  *   ini ;

    fXmolName=malloc(500*sizeof(char)); if (fXmolName==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): fXmolName[] malloc out of memory\n");   Error_no_free(errMsg); }
    fBoxSizeName=malloc(500*sizeof(char)); if (fBoxSizeName==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): fBoxSizeName[] malloc out of memory\n");   Error_no_free(errMsg); }

    ini = iniparser_load(filename);
    if (ini==NULL) {
        sprintf(errMsg,"Setup_ReadIniFile(): Error opening file %s",filename);  // Always test file open
        Error_no_free(errMsg);
    }
    
    //box
    ISNOTCUBIC = iniparser_getint(ini, "box:box_type", -1);
    strcpy(fBoxSizeName, (char*)iniparser_getstring(ini, "box:box_name", "-1"));    
    
    //run
    strcpy(fXmolName, (char*)iniparser_getstring(ini, "run:xyzfilename", "-1"));
    FRAMES = iniparser_getint(ini, "run:frames", -1);
    TOTALFRAMES=iniparser_getint(ini, "run:totalframes", -1);
    N = iniparser_getint(ini, "run:num_particles", -1);
    NA = iniparser_getint(ini, "run:numA_particles", -1);
    RHO = iniparser_getdouble(ini, "run:number_density", -1);
    TSTART = iniparser_getdouble(ini, "run:simulationstarttime", -1);
    FRAMETSTEP = iniparser_getdouble(ini, "run:simulationtimestep", -1);
    TFINAL = iniparser_getdouble(ini, "run:simulationendtime", -1);
    STARTFROM = iniparser_getint(ini, "run:start_from", -1);
    SAMPLEFREQ = iniparser_getint(ini, "run:sample_freqency", -1);

    //simulation
    rcutAA = iniparser_getdouble(ini, "simulation:rcutAA", -1);
    rcutAB = iniparser_getdouble(ini, "simulation:rcutAB", -1);
    rcutBB = iniparser_getdouble(ini, "simulation:rcutBB", -1);
    Vor = iniparser_getboolean(ini, "simulation:bond_type", -1);
    PBCs = iniparser_getboolean(ini, "simulation:pbcs", -1);
    fc = iniparser_getdouble(ini, "simulation:voronoi_parameter", -1);
    nB = iniparser_getint(ini, "simulation:num_bonds", -1);
    USELIST = iniparser_getboolean(ini, "simulation:cell_list", -1);

    //output
    doWriteBonds = iniparser_getboolean(ini, "output:bonds", -1);
    doWriteClus = iniparser_getboolean(ini, "output:clusts", -1);
    doWriteRaw = iniparser_getboolean(ini, "output:raw", -1);
    do11AcenXmol = iniparser_getboolean(ini, "output:11a", -1);
    do13AcenXmol = iniparser_getboolean(ini, "output:13a", -1);
    doWritePopPerFrame = iniparser_getboolean(ini, "output:pop_per_frame", -1);
    doSubClusts = iniparser_getboolean(ini, "output:subclusters", -1);
    PRINTINFO = iniparser_getboolean(ini, "extra:debug", -1);
    iniparser_getdouble(ini, "extra:shear", -1);

    // calculate derived values
    rcutAA2=rcutAA*rcutAA;
    rcutAB2=rcutAB*rcutAB;
    rcutBB2=rcutBB*rcutBB;
    if (Vor==1) {   // if using modified Voronoi method can't have different cut off lengths for the bonds
        rcutAB=rcutAA;
        rcutBB=rcutAA;
        rcutAB2=rcutAA2;
        rcutBB2=rcutAA2;
        printf("As voronoi method no individually specie-specie interaction length\nrcut %lg rcut2 %lg\n",rcutAA,rcutAA2);
    }
    initNoStatic=incrStatic=initNoClustPerPart=incrClustPerPart=1;

    if (NA>N) {
        Error_no_free("Setup_ReadIniFile(): NA > N - something wrong in xmol trajectory params file\n");
    }

    sidex=sidey=sidez=pow((double)N/RHO, 1.0/3.0);
    halfSidex=halfSidey=halfSidez=sidex/2.0;
    
    // print out values read from ini file
    printf("Xmol file name:%s Box file name:%s\n", fXmolName, fBoxSizeName);
    printf("ISNOTCUBIC %d\n",ISNOTCUBIC);
    printf("FRAMES %d N %d NA %d RHO %lg TSTART %lg\n",FRAMES,N,NA,RHO,TSTART);
    printf("FRAMETSTEP %lg TFINAL %lg STARTFROM %d SAMPLEFREQ %d\n",FRAMETSTEP,TFINAL, STARTFROM,SAMPLEFREQ);
    printf("rcutAA %lg rcutAB %lg rcutBB %lg\n",rcutAA,rcutAB,rcutBB);
    printf("rcutAA2 %lg rcutAB2 %lg rcutBB2 %lg\n",rcutAA2,rcutAB2,rcutBB2);
    printf("Vor %d PBCs %d fc %lg nB %d USELIST %d\n",Vor,PBCs,fc,nB,USELIST);
    printf("write bonds file %d doWriteClus %d doWriteRaw %d doWritePopPerFrame %d\n",doWriteBonds,doWriteClus,doWriteRaw,doWritePopPerFrame);
    printf("doSubClusts %d PRINTINFO %d\n",doSubClusts, PRINTINFO);

    if (ISNOTCUBIC==0) {
        printf("calculating box sides from RHO\n");
        printf("box side length = %.5lg, half side: %.5lg\n\n",sidex,halfSidex);
    }

    iniparser_freedict(ini);
}

void Setup_ReadBox(FILE *readIn)  {
    int sweep;
        
    if (feof(readIn)) Error("Setup_ReadBox(): end of input file reached\n");
    if(ISNOTCUBIC==3){
        printf("Triclinic Boundary Conditions\n");
      
        fscanf(readIn,"%i %lf %lf %lf %lf %lf %lf\n", &sweep,&sidex,&sidey,&sidez, &tiltxy, &tiltxz, &tiltyz);
        printf("iter Lx Ly Lz xy xz yz\n");

        printf("%i %lf %lf %lf %lf %lf %lf\n",sweep,sidex,sidey,sidez, tiltxy, tiltxz, tiltyz);}
    else
    {
        fscanf(readIn,"%i %lg %lg %lg\n", &sweep,&sidex,&sidey,&sidez);
        printf("%i %lg %lg %lg\n",sweep,sidex,sidey,sidez);
    }
    halfSidex = sidex/2.0;
    halfSidey = sidey/2.0;
    halfSidez = sidez/2.0;
}

void Setup_Readxyz(int e, int write, int f, FILE *readin) {     // read configuration from xmol trajectory
    int i;
    char c;
    double tx, ty, tz;
    int cntA;
    char input[1000],errMsg[1000];
    
    cntA=0; // check number of A-species is as expected
    
    if (feof(readin)) Error("Setup_Readxyz(): end of input file reached\n");
    fscanf(readin,"%d\n", &i);  // read number of particles
    if (i!=N) { // check number of particles is as expecting
        sprintf(errMsg,"Setup_Readxyz(): N %d from input frame %d does not match N %d from params file\n",i,e,N);
        Error(errMsg);
    }
    if (feof(readin)) Error("Setup_Readxyz(): end of input file reached\n");
    fgets(input,1000,readin);   // get the information line and do nothing with it
    if (write==1) { // if storing configurations, do some stuff
        if (PRINTINFO==1) printf("Setup_Readxyz(): TCC analysis frame %d/%d - reading in %d particles frame %d of XMOL\n",f,FRAMES-1,i,e);
    }
    for (i=0; i<N; ++i) {   // loop over all particles
        if (feof(readin)) Error("Setup_Readxyz(): end of input file reached\n");    // check not reached end of file
        fscanf(readin," %c  %lg %lg %lg\n", &c,&tx,&ty,&tz);
        if (c=='A') cntA++;     // check number of A-species is as expected
        else if (c=='C') cntA++;    // check number of A-species is as expected
        if (write==1) { // keeping configuration
            if (c=='A') rtype[i]=1; // set rtype array denoting cluster species
            else if (c=='B') rtype[i]=2;
            else if (c=='C') rtype[i]=1;
            else {
                sprintf(errMsg,"Setup_Readxyz(): unrecognized character of particle i %d from input frame %d\n",i,e);
                Error(errMsg);
            }
            if (PBCs==1 && ISNOTCUBIC!=3) {  // wrap particles back into the box
                if (tx<-halfSidex) { tx+=sidex; }
                else if (tx>halfSidex)   { tx-=sidex; }
                if (ty<-halfSidey) { ty+=sidey; }
                else if (ty>halfSidey)   { ty-=sidey; }
                if (tz<-halfSidez) { tz+=sidez; }
                else if (tz>halfSidez)   { tz-=sidez; }
            }
            x[i]=tx;    y[i]=ty;    z[i]=tz;    // set positions
            if (PRINTINFO==1) if (i==N-1) printf("f%d part%d %c %.5lg %.5lg %.5lg\n\n",f,i,c,x[i],y[i],z[i]);
        }
    }
    if (cntA!=NA) { // check number of A-species is as expected
        sprintf(errMsg,"Setup_Readxyz(): NA %d from input frame %d does not match NA %d from params file\n",cntA,e,NA);
        Error(errMsg);
    }
}

void Setup_InitStaticVars() { // Initialize lots of important variables for static TCC algorithm
    int f, j, k;
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
    
    if (USELIST==1) {   // if using cell list to detect bond network allocate memory
        head=malloc((ncells+1)*sizeof(int));    if (head==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): head[] malloc out of memory\n");  Error_no_free(errMsg); }
        map=malloc((13*ncells+1)*sizeof(int));  if (map==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): map[] malloc out of memory\n");    Error_no_free(errMsg); }
        llist=malloc((N+1)*sizeof(int));    if (llist==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): llist[] malloc out of memory\n");    Error_no_free(errMsg); }
        
        for (j=0; j<ncells+1; j++) head[j]=0;
        for (j=0; j<13*ncells+1; j++) map[j]=0;
        for (j=0; j<N+1; j++) llist[j]=0;
    }
    
    mean_pop_per_frame_sp3=mean_pop_per_frame_sp3a=mean_pop_per_frame_sp3b=mean_pop_per_frame_sp3c=0.0; // mean particle fractions of clusters
    mean_pop_per_frame_sp4=mean_pop_per_frame_sp4a=mean_pop_per_frame_sp4b=mean_pop_per_frame_sp4c=0.0;
    mean_pop_per_frame_sp5=mean_pop_per_frame_sp5a=mean_pop_per_frame_sp5b=mean_pop_per_frame_sp5c=0.0;
    mean_pop_per_frame_6Z=mean_pop_per_frame_7K=0.0;
    mean_pop_per_frame_8A=mean_pop_per_frame_8B=mean_pop_per_frame_8K=0.0;
    mean_pop_per_frame_9A=mean_pop_per_frame_9B=mean_pop_per_frame_9K=0.0;
    mean_pop_per_frame_10A=mean_pop_per_frame_10B=mean_pop_per_frame_10K=mean_pop_per_frame_10W=0.0;
    mean_pop_per_frame_11A=mean_pop_per_frame_11B=mean_pop_per_frame_11C=mean_pop_per_frame_11E=mean_pop_per_frame_11F=mean_pop_per_frame_11W=0.0;
    mean_pop_per_frame_12A=mean_pop_per_frame_12B=mean_pop_per_frame_12D=mean_pop_per_frame_12E=mean_pop_per_frame_12K=0.0;
    mean_pop_per_frame_13A=mean_pop_per_frame_13B=mean_pop_per_frame_13K=0.0;
    mean_pop_per_frame_FCC=mean_pop_per_frame_HCP=mean_pop_per_frame_BCC_9=mean_pop_per_frame_BCC_15=0.0;

    x = malloc(N*sizeof(double));   if (x==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): x[] malloc out of memory\n");    Error_no_free(errMsg); }    // positions of particles in a configuration
    y = malloc(N*sizeof(double));   if (y==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): y[] malloc out of memory\n");    Error_no_free(errMsg); }
    z = malloc(N*sizeof(double));   if (z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): z[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    rtype=malloc(N*sizeof(int)); if (rtype==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): rtype[] malloc out of memory\n");   Error_no_free(errMsg); }    // type of species

    raw_file_pointers = malloc(sizeof(FILE*) * num_cluster_types);
    cluster_file_pointers = malloc(sizeof(FILE*) * num_cluster_types);

    cnb = malloc(N*sizeof(int));    if (cnb==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): cnb[] malloc out of memory\n");    Error_no_free(errMsg); }    // number of "bonded" neighbours of a particle
    correctedBonds=0;   // count number of times have make j bonded to i given i bonded to j due to round off error in Voronoi code
    bNums = malloc(N*sizeof(int *));    if (bNums==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): bNums[] malloc out of memory\n");    Error_no_free(errMsg); }    // list of bonded particles to each particle
    for (j=0; j<N; ++j) { bNums[j] = malloc(nB*sizeof(int));    if (bNums[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): bNums[][] malloc out of memory\n");   Error_no_free(errMsg); } }

    bondlengths = malloc(N*sizeof(double *));   if (bondlengths==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): bondlengths[] malloc out of memory\n");    Error_no_free(errMsg); }    // array of bond lengths
    for (j=0; j<N; ++j) { bondlengths[j] = malloc(nB*sizeof(double));   if (bondlengths[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): bondlengths[][] malloc out of memory\n");   Error_no_free(errMsg); } }  
    
    for (j=0; j<N; j++) {   // reset the variables
        x[j]=y[j]=z[j]=0.0;
        rtype[j]=0;
        cnb[j]=0;
        
        for (k=0; k<nB; ++k) {
            bNums[j][k]=0;
            bondlengths[j][k]=0.0;
        }
    }

    maxnb=0;
    maxto3=maxto4=maxto5=0;
    
    // count number of bonds between spindle particles
    nsp3c_spindlebonds = malloc(FRAMES*sizeof(int));    if (nsp3c_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3c_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp4c_spindlebonds = malloc(FRAMES*sizeof(int));    if (nsp4c_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4c_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp5c_spindlebonds = malloc(FRAMES*sizeof(int));    if (nsp5c_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5c_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    // count number of times more than two spindles found to SP3/4/5 ring
    nsp3_excess_spindles = malloc(FRAMES*sizeof(int));  if (nsp3_excess_spindles==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3_excess_spindles[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp4_excess_spindles = malloc(FRAMES*sizeof(int));  if (nsp4_excess_spindles==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4_excess_spindles[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp5_excess_spindles = malloc(FRAMES*sizeof(int));  if (nsp5_excess_spindles==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5_excess_spindles[] malloc out of memory\n");  Error_no_free(errMsg); }
        
    // number of times each cluster type is found in each frame
    nsp3 = malloc(FRAMES*sizeof(int));  if (nsp3==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp3a = malloc(FRAMES*sizeof(int)); if (nsp3a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3a[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp3b = malloc(FRAMES*sizeof(int)); if (nsp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3b[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp3c = malloc(FRAMES*sizeof(int)); if (nsp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3c[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    nsp4 = malloc(FRAMES*sizeof(int));  if (nsp4==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp4a = malloc(FRAMES*sizeof(int)); if (nsp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4a[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp4b = malloc(FRAMES*sizeof(int)); if (nsp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4b[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp4c = malloc(FRAMES*sizeof(int)); if (nsp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4c[] malloc out of memory\n");    Error_no_free(errMsg); }

    nsp5 = malloc(FRAMES*sizeof(int));  if (nsp5==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp5a = malloc(FRAMES*sizeof(int)); if (nsp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5a[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp5b = malloc(FRAMES*sizeof(int)); if (nsp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5b[] malloc out of memory\n");    Error_no_free(errMsg); }
    nsp5c = malloc(FRAMES*sizeof(int)); if (nsp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp5c[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    n6Z = malloc(FRAMES*sizeof(int));   if (n6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n6Z[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    n7K = malloc(FRAMES*sizeof(int));   if (n7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n7K[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    n8A = malloc(FRAMES*sizeof(int));   if (n8A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n8A[] malloc out of memory\n");    Error_no_free(errMsg); }
    n8B = malloc(FRAMES*sizeof(int));   if (n8B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n8B[] malloc out of memory\n");    Error_no_free(errMsg); }
    n8K = malloc(FRAMES*sizeof(int));   if (n8K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n8K[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    n9A = malloc(FRAMES*sizeof(int));   if (n9A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n9A[] malloc out of memory\n");    Error_no_free(errMsg); }
    n9B = malloc(FRAMES*sizeof(int));   if (n9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n9B[] malloc out of memory\n");    Error_no_free(errMsg); }
    n9K = malloc(FRAMES*sizeof(int));   if (n9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n9K[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    n10A = malloc(FRAMES*sizeof(int));  if (n10A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n10A[] malloc out of memory\n");  Error_no_free(errMsg); }
    n10B = malloc(FRAMES*sizeof(int));  if (n10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n10B[] malloc out of memory\n");  Error_no_free(errMsg); }
    n10K = malloc(FRAMES*sizeof(int));  if (n10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n10K[] malloc out of memory\n");  Error_no_free(errMsg); }
    n10W = malloc(FRAMES*sizeof(int));  if (n10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n10W[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    n11A = malloc(FRAMES*sizeof(int));  if (n11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11A[] malloc out of memory\n");  Error_no_free(errMsg); }
    n11B = malloc(FRAMES*sizeof(int));  if (n11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11B[] malloc out of memory\n");  Error_no_free(errMsg); }
    n11C = malloc(FRAMES*sizeof(int));  if (n11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11C[] malloc out of memory\n");  Error_no_free(errMsg); }
    n11E = malloc(FRAMES*sizeof(int));  if (n11E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11E[] malloc out of memory\n");  Error_no_free(errMsg); }
    n11F = malloc(FRAMES*sizeof(int));  if (n11F==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11F[] malloc out of memory\n");  Error_no_free(errMsg); }
    n11W = malloc(FRAMES*sizeof(int));  if (n11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n11W[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    n12A = malloc(FRAMES*sizeof(int));  if (n12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n12A[] malloc out of memory\n");  Error_no_free(errMsg); }
    n12B = malloc(FRAMES*sizeof(int));  if (n12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n12B[] malloc out of memory\n");  Error_no_free(errMsg); }
    n12D = malloc(FRAMES*sizeof(int));  if (n12D==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n12D[] malloc out of memory\n");  Error_no_free(errMsg); }
    n12E = malloc(FRAMES*sizeof(int));  if (n12E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n12E[] malloc out of memory\n");  Error_no_free(errMsg); }
    n12K = malloc(FRAMES*sizeof(int));  if (n12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n12K[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    n13A = malloc(FRAMES*sizeof(int));  if (n13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n13A[] malloc out of memory\n");  Error_no_free(errMsg); }
    n13B = malloc(FRAMES*sizeof(int));  if (n13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n13B[] malloc out of memory\n");  Error_no_free(errMsg); }
    n13K = malloc(FRAMES*sizeof(int));  if (n13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n13K[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    nFCC = malloc(FRAMES*sizeof(int));  if (nFCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nFCC[] malloc out of memory\n");  Error_no_free(errMsg); }
    nHCP = malloc(FRAMES*sizeof(int));  if (nHCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nHCP[] malloc out of memory\n");  Error_no_free(errMsg); }
    nBCC_9 = malloc(FRAMES*sizeof(int));    if (nBCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nBCC_9[] malloc out of memory\n");  Error_no_free(errMsg); }
    nBCC_15 = malloc(FRAMES*sizeof(int));   if (nBCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nBCC_15[] malloc out of memory\n");    Error_no_free(errMsg); }
    
    for (f=0; f<FRAMES; f++) { // reset the variables
        nsp3c_spindlebonds[f]=0;
        nsp4c_spindlebonds[f]=0;
        nsp5c_spindlebonds[f]=0;
        
        nsp3_excess_spindles[f]=0;
        nsp4_excess_spindles[f]=0;
        nsp5_excess_spindles[f]=0;
        
        nsp3[f]=0;
        nsp3a[f]=0;
        nsp3b[f]=0;
        nsp3c[f]=0;
        
        nsp4[f]=0;
        nsp4a[f]=0;
        nsp4b[f]=0;
        nsp4c[f]=0;

        nsp5[f]=0;
        nsp5a[f]=0;
        nsp5b[f]=0;
        nsp5c[f]=0;
        
        n6Z[f]=0;
        
        n7K[f]=0;
                
        n8A[f]=0;
        n8B[f]=0;
        n8K[f]=0;
        
        n9A[f]=0;
        n9B[f]=0;
        n9K[f]=0;
        
        n10A[f]=0;
        n10B[f]=0;
        n10K[f]=0;
        n10W[f]=0;
        
        n11A[f]=0;
        n11B[f]=0;
        n11C[f]=0;
        n11E[f]=0;
        n11F[f]=0;
        n11W[f]=0;
        
        n12A[f]=0;
        n12B[f]=0;
        n12D[f]=0;
        n12E[f]=0;
        n12K[f]=0;
        
        n13A[f]=0;
        n13B[f]=0;
        n13K[f]=0;
        
        nFCC[f]=0;
        nHCP[f]=0;
        nBCC_9[f]=0;
        nBCC_15[f]=0;
    }
    
    // arrays for the clusters found in each frame
    sp3a = malloc(msp3a*sizeof(int *)); if (sp3a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3a; ++j) { sp3a[j] = malloc(3*sizeof(int));  if (sp3a[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp3b = malloc(msp3b*sizeof(int *)); if (sp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3b; ++j) { sp3b[j] = malloc(4*sizeof(int));  if (sp3b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp3c = malloc(msp3c*sizeof(int *)); if (sp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp3c; ++j) { sp3c[j] = malloc(5*sizeof(int));  if (sp3c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp3c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp4a = malloc(msp4a*sizeof(int *)); if (sp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4a; ++j) { sp4a[j] = malloc(4*sizeof(int));  if (sp4a[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp4b = malloc(msp4b*sizeof(int *)); if (sp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4b; ++j) { sp4b[j] = malloc(5*sizeof(int));  if (sp4b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp4c = malloc(msp4c*sizeof(int *)); if (sp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4c; ++j) { sp4c[j] = malloc(6*sizeof(int));  if (sp4c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp4c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp5a = malloc(msp5a*sizeof(int *)); if (sp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5a[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5a; ++j) { sp5a[j] = malloc(5*sizeof(int));  if (sp5a[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5a[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp5b = malloc(msp5b*sizeof(int *)); if (sp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5b; ++j) { sp5b[j] = malloc(6*sizeof(int));  if (sp5b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    sp5c = malloc(msp5c*sizeof(int *)); if (sp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp5c; ++j) { sp5c[j] = malloc(7*sizeof(int));  if (sp5c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sp5c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    // arrays for the number of clusters of each type bonded to each particle
    mem_sp3b = malloc(N*sizeof(int *)); if (mem_sp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp3b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp3b[j] = malloc(mmem_sp3b*sizeof(int));  if (mem_sp3b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp3b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp3b = malloc(N*sizeof(int));  if (nmem_sp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp3b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp3c = malloc(N*sizeof(int *)); if (mem_sp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp3c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp3c[j] = malloc(mmem_sp3c*sizeof(int));  if (mem_sp3c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp3c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp3c = malloc(N*sizeof(int));  if (nmem_sp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp3c[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp4b = malloc(N*sizeof(int *)); if (mem_sp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp4b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp4b[j] = malloc(mmem_sp4b*sizeof(int));  if (mem_sp4b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp4b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp4b = malloc(N*sizeof(int));  if (nmem_sp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp4b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp4c = malloc(N*sizeof(int *)); if (mem_sp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp4c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp4c[j] = malloc(mmem_sp4c*sizeof(int));  if (mem_sp4c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp4c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp4c = malloc(N*sizeof(int));  if (nmem_sp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp4c[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp5b = malloc(N*sizeof(int *)); if (mem_sp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp5b[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp5b[j] = malloc(mmem_sp5b*sizeof(int));  if (mem_sp5b[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp5b[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp5b = malloc(N*sizeof(int));  if (nmem_sp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp5b[] malloc out of memory\n");    Error_no_free(errMsg); }
    mem_sp5c = malloc(N*sizeof(int *)); if (mem_sp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp5c[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<N; ++j) { mem_sp5c[j] = malloc(mmem_sp5c*sizeof(int));  if (mem_sp5c[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): mem_sp5c[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    nmem_sp5c = malloc(N*sizeof(int));  if (nmem_sp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nmem_sp5c[] malloc out of memory\n");    Error_no_free(errMsg); }

    // arrays for the clusters found in each frame
    hc6Z = malloc(m6Z*sizeof(int *));   if (hc6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc6Z[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m6Z; ++j) { hc6Z[j] = malloc(6*sizeof(int));    if (hc6Z[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc6Z[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc7K = malloc(m7K*sizeof(int *));   if (hc7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc7K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m7K; ++j) { hc7K[j] = malloc(7*sizeof(int));    if (hc7K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc7K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc8A = malloc(m8A*sizeof(int *));   if (hc8A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8A[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8A; ++j) { hc8A[j] = malloc(8*sizeof(int));    if (hc8A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8A[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc8B = malloc(m8B*sizeof(int *));   if (hc8B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8B[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8B; ++j) { hc8B[j] = malloc(8*sizeof(int));    if (hc8B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8B[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc8K = malloc(m8K*sizeof(int *));   if (hc8K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m8K; ++j) { hc8K[j] = malloc(8*sizeof(int));    if (hc8K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc8K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc9A = malloc(m9A*sizeof(int *));   if (hc9A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9A[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9A; ++j) { hc9A[j] = malloc(9*sizeof(int));    if (hc9A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9A[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc9B = malloc(m9B*sizeof(int *));   if (hc9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9B[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9B; ++j) { hc9B[j] = malloc(9*sizeof(int));    if (hc9B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9B[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    hc9K = malloc(m9K*sizeof(int *));   if (hc9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9K[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<m9K; ++j) { hc9K[j] = malloc(9*sizeof(int));    if (hc9K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc9K[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    hc10A = malloc(m10A*sizeof(int *)); if (hc10A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10A; ++j) { hc10A[j] = malloc(10*sizeof(int)); if (hc10A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10B = malloc(m10B*sizeof(int *)); if (hc10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10B; ++j) { hc10B[j] = malloc(10*sizeof(int)); if (hc10B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10K = malloc(m10K*sizeof(int *)); if (hc10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10K; ++j) { hc10K[j] = malloc(10*sizeof(int)); if (hc10K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc10W = malloc(m10W*sizeof(int *)); if (hc10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10W[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m10W; ++j) { hc10W[j] = malloc(10*sizeof(int)); if (hc10W[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc10W[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hc11A = malloc(m11A*sizeof(int *)); if (hc11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11A; ++j) { hc11A[j] = malloc(11*sizeof(int)); if (hc11A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11B = malloc(m11B*sizeof(int *)); if (hc11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11B; ++j) { hc11B[j] = malloc(11*sizeof(int)); if (hc11B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11C = malloc(m11C*sizeof(int *)); if (hc11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11C[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11C; ++j) { hc11C[j] = malloc(11*sizeof(int)); if (hc11C[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11C[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11E = malloc(m11E*sizeof(int *)); if (hc11E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11E[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11E; ++j) { hc11E[j] = malloc(11*sizeof(int)); if (hc11E[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11E[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11F = malloc(m11F*sizeof(int *)); if (hc11F==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11F[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11F; ++j) { hc11F[j] = malloc(11*sizeof(int)); if (hc11F[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11F[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc11W = malloc(m11W*sizeof(int *)); if (hc11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11W[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m11W; ++j) { hc11W[j] = malloc(11*sizeof(int)); if (hc11W[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc11W[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hc12A = malloc(m12A*sizeof(int *)); if (hc12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12A; ++j) { hc12A[j] = malloc(12*sizeof(int)); if (hc12A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12B = malloc(m12B*sizeof(int *)); if (hc12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12B; ++j) { hc12B[j] = malloc(12*sizeof(int)); if (hc12B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12D = malloc(m12D*sizeof(int *)); if (hc12D==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12D[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12D; ++j) { hc12D[j] = malloc(12*sizeof(int)); if (hc12D[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12D[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12E = malloc(m12E*sizeof(int *)); if (hc12E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12E[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12E; ++j) { hc12E[j] = malloc(12*sizeof(int)); if (hc12E[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12E[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc12K = malloc(m12K*sizeof(int *)); if (hc12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m12K; ++j) { hc12K[j] = malloc(12*sizeof(int)); if (hc12K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc12K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
        
    hc13A = malloc(m13A*sizeof(int *)); if (hc13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13A[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13A; ++j) { hc13A[j] = malloc(13*sizeof(int)); if (hc13A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13A[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc13B = malloc(m13B*sizeof(int *)); if (hc13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13B[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13B; ++j) { hc13B[j] = malloc(13*sizeof(int)); if (hc13B[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13B[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hc13K = malloc(m13K*sizeof(int *)); if (hc13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13K[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<m13K; ++j) { hc13K[j] = malloc(13*sizeof(int)); if (hc13K[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc13K[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    
    hcFCC = malloc(mFCC*sizeof(int *)); if (hcFCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcFCC[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mFCC; ++j) { hcFCC[j] = malloc(13*sizeof(int)); if (hcFCC[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcFCC[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcHCP = malloc(mHCP*sizeof(int *)); if (hcHCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcHCP[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mHCP; ++j) { hcHCP[j] = malloc(13*sizeof(int)); if (hcHCP[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcHCP[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcBCC_9 = malloc(mBCC_9*sizeof(int *)); if (hcBCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcBCC_9[] malloc out of memory\n");    Error_no_free(errMsg); }
    for (j=0; j<mBCC_9; ++j) { hcBCC_9[j] = malloc(9*sizeof(int));  if (hcBCC_9[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcBCC_9[][] malloc out of memory\n");   Error_no_free(errMsg); } }
    hcBCC_15 = malloc(mBCC_15*sizeof(int *));   if (hcBCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcBCC_15[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<mBCC_15; ++j) { hcBCC_15[j] = malloc(15*sizeof(int));   if (hcBCC_15[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hcBCC_15[][] malloc out of memory\n"); Error_no_free(errMsg); } }
    
    // reset the variables
    for (j=0; j<msp3a; ++j) { for (k=0; k<3; k++) { sp3a[j][k]=-1; }}
    for (j=0; j<msp3b; ++j) { for (k=0; k<4; k++) { sp3b[j][k]=-1; }}
    for (j=0; j<msp3c; ++j) { for (k=0; k<5; k++) { sp3c[j][k]=-1; }}
    
    for (j=0; j<N; ++j) { 
        for (k=0; k<mmem_sp3b; k++) mem_sp3b[j][k]=-1; 
        for (k=0; k<mmem_sp3c; k++) mem_sp3c[j][k]=-1; 
        for (k=0; k<mmem_sp4b; k++) mem_sp4b[j][k]=-1; 
        for (k=0; k<mmem_sp4c; k++) mem_sp4c[j][k]=-1;
        for (k=0; k<mmem_sp5b; k++) mem_sp5b[j][k]=-1; 
        for (k=0; k<mmem_sp5c; k++) mem_sp5c[j][k]=-1;
        nmem_sp3b[j]=0;
        nmem_sp3c[j]=0;
        nmem_sp4b[j]=0;
        nmem_sp4c[j]=0;
        nmem_sp5b[j]=0;
        nmem_sp5c[j]=0;
    }
    
    for (j=0; j<msp4a; ++j) { for (k=0; k<4; k++) { sp4a[j][k]=-1; }}
    for (j=0; j<msp4b; ++j) { for (k=0; k<5; k++) { sp4b[j][k]=-1; }}
    for (j=0; j<msp4c; ++j) { for (k=0; k<6; k++) { sp4c[j][k]=-1; }}
        
    for (j=0; j<msp5a; ++j) { for (k=0; k<5; k++) { sp5a[j][k]=-1; }}
    for (j=0; j<msp5b; ++j) { for (k=0; k<6; k++) { sp5b[j][k]=-1; }}
    for (j=0; j<msp5c; ++j) { for (k=0; k<7; k++) { sp5c[j][k]=-1; }}
    
    for (j=0; j<m6Z; ++j) { for (k=0; k<6; k++) { hc6Z[j][k]=-1; }}
    
    for (j=0; j<m7K; ++j) { for (k=0; k<7; k++) { hc7K[j][k]=-1; }}
    
    for (j=0; j<m8A; ++j) { for (k=0; k<8; k++) { hc8A[j][k]=-1; }}
    for (j=0; j<m8B; ++j) { for (k=0; k<8; k++) { hc8B[j][k]=-1; }}
    for (j=0; j<m8K; ++j) { for (k=0; k<8; k++) { hc8K[j][k]=-1; }}
    
    for (j=0; j<m9A; ++j) { for (k=0; k<9; k++) { hc9A[j][k]=-1; }}
    for (j=0; j<m9B; ++j) { for (k=0; k<9; k++) { hc9B[j][k]=-1; }}
    for (j=0; j<m9K; ++j) { for (k=0; k<9; k++) { hc9K[j][k]=-1; }}
    
    for (j=0; j<m10A; ++j) { for (k=0; k<10; k++) { hc10A[j][k]=-1; }}
    for (j=0; j<m10B; ++j) { for (k=0; k<10; k++) { hc10B[j][k]=-1; }}
    for (j=0; j<m10K; ++j) { for (k=0; k<10; k++) { hc10K[j][k]=-1; }}
    for (j=0; j<m10W; ++j) { for (k=0; k<10; k++) { hc10W[j][k]=-1; }}
    
    for (j=0; j<m11A; ++j) { for (k=0; k<11; k++) { hc11A[j][k]=-1; }}
    for (j=0; j<m11B; ++j) { for (k=0; k<11; k++) { hc11B[j][k]=-1; }}
    for (j=0; j<m11C; ++j) { for (k=0; k<11; k++) { hc11C[j][k]=-1; }}
    for (j=0; j<m11E; ++j) { for (k=0; k<11; k++) { hc11E[j][k]=-1; }}
    for (j=0; j<m11F; ++j) { for (k=0; k<11; k++) { hc11F[j][k]=-1; }}
    for (j=0; j<m11W; ++j) { for (k=0; k<11; k++) { hc11W[j][k]=-1; }}
    
    for (j=0; j<m12A; ++j) { for (k=0; k<12; k++) { hc12A[j][k]=-1; }}
    for (j=0; j<m12B; ++j) { for (k=0; k<12; k++) { hc12B[j][k]=-1; }}
    for (j=0; j<m12D; ++j) { for (k=0; k<12; k++) { hc12D[j][k]=-1; }}
    for (j=0; j<m12E; ++j) { for (k=0; k<12; k++) { hc12E[j][k]=-1; }}
    for (j=0; j<m12K; ++j) { for (k=0; k<12; k++) { hc12K[j][k]=-1; }}
    
    for (j=0; j<m13A; ++j) { for (k=0; k<13; k++) { hc13A[j][k]=-1; }}
    for (j=0; j<m13B; ++j) { for (k=0; k<13; k++) { hc13B[j][k]=-1; }}
    for (j=0; j<m13K; ++j) { for (k=0; k<13; k++) { hc13K[j][k]=-1; }}
    
    for (j=0; j<mFCC; ++j) { for (k=0; k<13; k++) { hcFCC[j][k]=-1; }}
    for (j=0; j<mHCP; ++j) { for (k=0; k<13; k++) { hcHCP[j][k]=-1; }}
    for (j=0; j<mBCC_9; ++j) { for (k=0; k<9; k++) { hcBCC_9[j][k]=-1; }}
    for (j=0; j<mBCC_15; ++j) { for (k=0; k<15; k++) { hcBCC_15[j][k]=-1; }}
    
    
    // character arrays listing what type each particle is when found in a cluster
    ssp3=malloc(N*sizeof(char)); if (ssp3==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp3[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp3a=malloc(N*sizeof(char)); if (ssp3a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp3a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp3b=malloc(N*sizeof(char)); if (ssp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp3b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp3c=malloc(N*sizeof(char)); if (ssp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp3c[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp4=malloc(N*sizeof(char)); if (ssp4==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4a=malloc(N*sizeof(char)); if (ssp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4b=malloc(N*sizeof(char)); if (ssp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4c=malloc(N*sizeof(char)); if (ssp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4c[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s6Z=malloc(N*sizeof(char)); if (s6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s6Z[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s7K=malloc(N*sizeof(char)); if (s7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s7K[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp5=malloc(N*sizeof(char)); if (ssp5==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5a=malloc(N*sizeof(char)); if (ssp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5b=malloc(N*sizeof(char)); if (ssp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5b[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5c=malloc(N*sizeof(char)); if (ssp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5c[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s8A=malloc(N*sizeof(char)); if (s8A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s8A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s8B=malloc(N*sizeof(char)); if (s8B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s8B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s8K=malloc(N*sizeof(char)); if (s8K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s8K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9A=malloc(N*sizeof(char)); if (s9A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9B=malloc(N*sizeof(char)); if (s9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9K=malloc(N*sizeof(char)); if (s9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10A=malloc(N*sizeof(char)); if (s10A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10B=malloc(N*sizeof(char)); if (s10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10K=malloc(N*sizeof(char)); if (s10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10W=malloc(N*sizeof(char)); if (s10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10W[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11A=malloc(N*sizeof(char)); if (s11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11B=malloc(N*sizeof(char)); if (s11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11C=malloc(N*sizeof(char)); if (s11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11C[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11E=malloc(N*sizeof(char)); if (s11E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11E[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11F=malloc(N*sizeof(char)); if (s11F==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11F[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11W=malloc(N*sizeof(char)); if (s11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11W[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12A=malloc(N*sizeof(char)); if (s12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12B=malloc(N*sizeof(char)); if (s12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12D=malloc(N*sizeof(char)); if (s12D==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12D[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12E=malloc(N*sizeof(char)); if (s12E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12E[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12K=malloc(N*sizeof(char)); if (s12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12K[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13A=malloc(N*sizeof(char)); if (s13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13A[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13B=malloc(N*sizeof(char)); if (s13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13B[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13K=malloc(N*sizeof(char)); if (s13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13K[] malloc out of memory\n"); Error_no_free(errMsg); }
    sFCC=malloc(N*sizeof(char)); if (sFCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sFCC[] malloc out of memory\n"); Error_no_free(errMsg); }
    sHCP=malloc(N*sizeof(char)); if (sHCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sHCP[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_9=malloc(N*sizeof(char)); if (sBCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_9[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_15=malloc(N*sizeof(char)); if (sBCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_15[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s9B_cen=malloc(N*sizeof(char)); if (s9B_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9B_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9K_cen=malloc(N*sizeof(char)); if (s9K_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9K_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10B_cen=malloc(N*sizeof(char)); if (s10B_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10B_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10K_cen=malloc(N*sizeof(char)); if (s10K_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10K_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10W_cen=malloc(N*sizeof(char)); if (s10W_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10W_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11A_cen=malloc(N*sizeof(char)); if (s11A_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11A_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11B_cen=malloc(N*sizeof(char)); if (s11B_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11B_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11C_cen=malloc(N*sizeof(char)); if (s11C_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11C_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11W_cen=malloc(N*sizeof(char)); if (s11W_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11W_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12A_cen=malloc(N*sizeof(char)); if (s12A_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12A_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12B_cen=malloc(N*sizeof(char)); if (s12B_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12B_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12K_cen=malloc(N*sizeof(char)); if (s12K_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12K_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13A_cen=malloc(N*sizeof(char)); if (s13A_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13A_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13B_cen=malloc(N*sizeof(char)); if (s13B_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13B_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13K_cen=malloc(N*sizeof(char)); if (s13K_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13K_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    sFCC_cen=malloc(N*sizeof(char)); if (sFCC_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sFCC_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    sHCP_cen=malloc(N*sizeof(char)); if (sHCP_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sHCP_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_9_cen=malloc(N*sizeof(char)); if (sBCC_9_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_9_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_15_cen=malloc(N*sizeof(char)); if (sBCC_15_cen==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_15_cen[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s9B_shell=malloc(N*sizeof(char)); if (s9B_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9B_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s9K_shell=malloc(N*sizeof(char)); if (s9K_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s9K_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10B_shell=malloc(N*sizeof(char)); if (s10B_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10B_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10K_shell=malloc(N*sizeof(char)); if (s10K_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10K_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s10W_shell=malloc(N*sizeof(char)); if (s10W_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s10W_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11A_shell=malloc(N*sizeof(char)); if (s11A_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11A_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11B_shell=malloc(N*sizeof(char)); if (s11B_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11B_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11C_shell=malloc(N*sizeof(char)); if (s11C_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11C_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s11W_shell=malloc(N*sizeof(char)); if (s11W_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s11W_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12A_shell=malloc(N*sizeof(char)); if (s12A_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12A_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12B_shell=malloc(N*sizeof(char)); if (s12B_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12B_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s12K_shell=malloc(N*sizeof(char)); if (s12K_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s12K_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13A_shell=malloc(N*sizeof(char)); if (s13A_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13A_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13B_shell=malloc(N*sizeof(char)); if (s13B_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13B_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    s13K_shell=malloc(N*sizeof(char)); if (s13K_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s13K_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    sFCC_shell=malloc(N*sizeof(char)); if (sFCC_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sFCC_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    sHCP_shell=malloc(N*sizeof(char)); if (sHCP_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sHCP_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_9_shell=malloc(N*sizeof(char)); if (sBCC_9_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_9_shell[] malloc out of memory\n"); Error_no_free(errMsg); }
    sBCC_15_shell=malloc(N*sizeof(char)); if (sBCC_15_shell==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): sBCC_15_shell[] malloc out of memory\n"); Error_no_free(errMsg); }

    for (j=0; j<N; j++) { // reset these variables
        ssp3[j]='C';
        ssp3a[j]='C';
        ssp3b[j]='C';
        ssp3c[j]='C';
        
        ssp4[j]='C';
        ssp4a[j]='C';
        ssp4b[j]='C';
        ssp4c[j]='C';
        
        ssp5[j]='C';
        ssp5a[j]='C';
        ssp5b[j]='C';
        ssp5c[j]='C';
        
        s6Z[j]='C';
        s7K[j]='C';
        s8A[j]='C';
        s8B[j]='C';
        s8K[j]='C';
        s9A[j]='C';
        s9B[j]='C';
        s9K[j]='C';
        s10A[j]='C';
        s10B[j]='C';
        s10K[j]='C';
        s10W[j]='C';
        s11A[j]='C';
        s11B[j]='C';
        s11C[j]='C';
        s11E[j]='C';
        s11F[j]='C';
        s11W[j]='C';
        s12A[j]='C';
        s12B[j]='C';
        s12D[j]='C';
        s12E[j]='C';
        s12K[j]='C';
        s13A[j]='C';
        s13B[j]='C';
        s13K[j]='C';
        sFCC[j]='C';
        sHCP[j]='C';
        sBCC_9[j]='C';
        sBCC_15[j]='C';
        
        s9B_cen[j]='C';
        s9K_cen[j]='C';
        s10B_cen[j]='C';
        s10K_cen[j]='C';
        s10W_cen[j]='C';
        s11A_cen[j]='C';
        s11B_cen[j]='C';
        s11C_cen[j]='C';
        s11W_cen[j]='C';
        s12A_cen[j]='C';
        s12B_cen[j]='C';
        s12K_cen[j]='C';
        s13A_cen[j]='C';
        s13B_cen[j]='C';
        s13K_cen[j]='C';
        sFCC_cen[j]='C';
        sHCP_cen[j]='C';
        sBCC_9_cen[j]='C';
        sBCC_15_cen[j]='C';
        
        s9B_shell[j]='C';
        s9K_shell[j]='C';
        s10B_shell[j]='C';
        s10K_shell[j]='C';
        s10W_shell[j]='C';
        s11A_shell[j]='C';
        s11B_shell[j]='C';
        s11C_shell[j]='C';
        s11W_shell[j]='C';
        s12A_shell[j]='C';
        s12B_shell[j]='C';
        s12K_shell[j]='C';
        s13A_shell[j]='C';
        s13B_shell[j]='C';
        s13K_shell[j]='C';
        sFCC_shell[j]='C';
        sHCP_shell[j]='C';
        sBCC_9_shell[j]='C';
        sBCC_15_shell[j]='C';
    }
    
    // particle fraction of particles in each cluster in each frame
    pop_per_frame_sp3 = malloc(FRAMES*sizeof(double));  if (pop_per_frame_sp3==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp3[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_sp3a = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp3a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp3a[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp3b = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp3b[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp3c = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp3c[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    pop_per_frame_sp4 = malloc(FRAMES*sizeof(double));  if (pop_per_frame_sp4==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp4[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_sp4a = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp4a[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp4b = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp4b[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp4c = malloc(FRAMES*sizeof(double));   if (pop_per_frame_sp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp4c[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    pop_per_frame_sp5 = malloc(FRAMES*sizeof(double));  if (pop_per_frame_sp5==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp5[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_sp5a = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp5a[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp5b = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp5b[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_sp5c = malloc(FRAMES*sizeof(double)); if (pop_per_frame_sp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_sp5c[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    pop_per_frame_6Z = malloc(FRAMES*sizeof(double));   if (pop_per_frame_6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_6Z[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_7K = malloc(FRAMES*sizeof(double));   if (pop_per_frame_7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_7K[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_8A = malloc(FRAMES*sizeof(double));   if (pop_per_frame_8A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_8A[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_8B = malloc(FRAMES*sizeof(double));   if (pop_per_frame_8B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_8B[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_8K = malloc(FRAMES*sizeof(double));   if (pop_per_frame_8K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_8K[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_9A = malloc(FRAMES*sizeof(double));   if (pop_per_frame_9A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_9A[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_9B = malloc(FRAMES*sizeof(double));   if (pop_per_frame_9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_9B[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_9K = malloc(FRAMES*sizeof(double));   if (pop_per_frame_9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_9K[] malloc out of memory\n");  Error_no_free(errMsg); }
    pop_per_frame_10A = malloc(FRAMES*sizeof(double));  if (pop_per_frame_10A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_10A[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_10B = malloc(FRAMES*sizeof(double));  if (pop_per_frame_10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_10B[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_10K = malloc(FRAMES*sizeof(double));  if (pop_per_frame_10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_10K[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_10W = malloc(FRAMES*sizeof(double));  if (pop_per_frame_10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_10W[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11A = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11A[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11B = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11B[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11C = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11C[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11E = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11E[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11F = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11F==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11F[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_11W = malloc(FRAMES*sizeof(double));  if (pop_per_frame_11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_11W[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_12A = malloc(FRAMES*sizeof(double));  if (pop_per_frame_12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_12A[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_12B = malloc(FRAMES*sizeof(double));  if (pop_per_frame_12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_12B[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_12D = malloc(FRAMES*sizeof(double));  if (pop_per_frame_12D==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_12D[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_12E = malloc(FRAMES*sizeof(double));  if (pop_per_frame_12E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_12E[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_12K = malloc(FRAMES*sizeof(double));  if (pop_per_frame_12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_12K[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_13A = malloc(FRAMES*sizeof(double));  if (pop_per_frame_13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_13A[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_13B = malloc(FRAMES*sizeof(double));  if (pop_per_frame_13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_13B[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_13K = malloc(FRAMES*sizeof(double));  if (pop_per_frame_13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_13K[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_FCC = malloc(FRAMES*sizeof(double));  if (pop_per_frame_FCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_FCC[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_HCP = malloc(FRAMES*sizeof(double));  if (pop_per_frame_HCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_HCP[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_BCC_9 = malloc(FRAMES*sizeof(double));    if (pop_per_frame_BCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_BCC_9[] malloc out of memory\n");    Error_no_free(errMsg); }
    pop_per_frame_BCC_15 = malloc(FRAMES*sizeof(double));   if (pop_per_frame_BCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_BCC_15[] malloc out of memory\n");  Error_no_free(errMsg); }
    
    for (f=0; f<FRAMES; f++) {  // reset these arrays
        pop_per_frame_sp3[f]=0;
        pop_per_frame_sp3a[f]=0;
        pop_per_frame_sp3b[f]=0;
        pop_per_frame_sp3c[f]=0;
        
        pop_per_frame_sp4[f]=0;
        pop_per_frame_sp4a[f]=0;
        pop_per_frame_sp4b[f]=0;
        pop_per_frame_sp4c[f]=0;
        
        pop_per_frame_sp5[f]=0;
        pop_per_frame_sp5a[f]=0;
        pop_per_frame_sp5b[f]=0;
        pop_per_frame_sp5c[f]=0;
        
        pop_per_frame_6Z[f]=0;
        pop_per_frame_7K[f]=0;
        pop_per_frame_8A[f]=0;
        pop_per_frame_8B[f]=0;
        pop_per_frame_8K[f]=0;
        pop_per_frame_9A[f]=0;
        pop_per_frame_9B[f]=0;
        pop_per_frame_9K[f]=0;
        pop_per_frame_10A[f]=0;
        pop_per_frame_10B[f]=0;
        pop_per_frame_10K[f]=0;
        pop_per_frame_10W[f]=0;
        pop_per_frame_11A[f]=0;
        pop_per_frame_11B[f]=0;
        pop_per_frame_11C[f]=0;
        pop_per_frame_11E[f]=0;
        pop_per_frame_11F[f]=0;
        pop_per_frame_11W[f]=0;
        pop_per_frame_12A[f]=0;
        pop_per_frame_12B[f]=0;
        pop_per_frame_12D[f]=0;
        pop_per_frame_12E[f]=0;
        pop_per_frame_12K[f]=0;
        pop_per_frame_13A[f]=0;
        pop_per_frame_13B[f]=0;
        pop_per_frame_13K[f]=0;
        pop_per_frame_FCC[f]=0;
        pop_per_frame_HCP[f]=0;
        pop_per_frame_BCC_9[f]=0;
        pop_per_frame_BCC_15[f]=0;
    }
}

void Setup_ResetStaticVars(int f) { // Reset static variables in each frame
    int i, j, k;

    for (i=0; i<N; ++i) {
        x[i]=y[i]=z[i]=-1000.0;
        for (j=0;j<nB;++j) {
            bNums[i][j]=-1;
            bondlengths[i][j]=-1.0;
        }
        cnb[i]=0;
    }
    
    nsp3c_spindlebonds[f]=nsp4c_spindlebonds[f]=nsp5c_spindlebonds[f]=0;
    nsp3_excess_spindles[f]=nsp4_excess_spindles[f]=nsp5_excess_spindles[f]=0;
    
    nsp3[f]=nsp3a[f]=nsp3b[f]=nsp3c[f]=0;
    nsp4[f]=nsp4a[f]=nsp4b[f]=nsp4c[f]=0;
    nsp5[f]=nsp5a[f]=nsp5b[f]=nsp5c[f]=0;
    
    n6Z[f]=n7K[f]=0;
    n8A[f]=n8B[f]=n8K[f]=0;
    n9A[f]=n9B[f]=n9K[f]=0;
    n10A[f]=n10B[f]=n10K[f]=n10W[f]=0;
    n11A[f]=n11B[f]=n11C[f]=n11E[f]=n11F[f]=n11W[f]=0;
    n12A[f]=n12B[f]=n12D[f]=n12E[f]=n12K[f]=0;
    n13A[f]=n13B[f]=n13K[f]=0;
    nFCC[f]=nHCP[f]=nBCC_9[f]=nBCC_15[f]=0;
    
    pop_per_frame_sp3[f]=pop_per_frame_sp3a[f]=pop_per_frame_sp3b[f]=pop_per_frame_sp3c[f]=0.0;
    pop_per_frame_sp4[f]=pop_per_frame_sp4a[f]=pop_per_frame_sp4b[f]=0.0;
    pop_per_frame_sp5[f]=pop_per_frame_sp5a[f]=pop_per_frame_sp5b[f]=pop_per_frame_sp5c[f]=0.0;
    pop_per_frame_sp4c[f]=pop_per_frame_6Z[f]=pop_per_frame_7K[f]=0.0;
    pop_per_frame_8A[f]=pop_per_frame_8B[f]=pop_per_frame_8K[f]=0.0;
    pop_per_frame_9A[f]=pop_per_frame_9B[f]=pop_per_frame_9K[f]=0.0;
    pop_per_frame_10A[f]=pop_per_frame_10B[f]=pop_per_frame_10K[f]=pop_per_frame_10W[f]=0.0;
    pop_per_frame_11A[f]=pop_per_frame_11B[f]=pop_per_frame_11C[f]=pop_per_frame_11E[f]=pop_per_frame_11F[f]=pop_per_frame_11W[f]=0.0;
    pop_per_frame_12A[f]=pop_per_frame_12B[f]=pop_per_frame_12D[f]=pop_per_frame_12E[f]=pop_per_frame_12K[f]=0.0;
    pop_per_frame_13A[f]=pop_per_frame_13B[f]=pop_per_frame_13K[f]=0.0;
    pop_per_frame_FCC[f]=pop_per_frame_HCP[f]=pop_per_frame_BCC_9[f]=pop_per_frame_BCC_15[f]=0.0;
    

    for (i=0; i<msp3a; ++i) {
        for (j=0;j<3;++j) sp3a[i][j]=-1;
    }
    for (i=0; i<msp3b; ++i) {
        for (j=0;j<4;++j) sp3b[i][j]=-1;
    }
    for (i=0; i<msp3c; ++i) {
        for (j=0;j<5;++j) sp3c[i][j]=-1;
    }
    for (j=0; j<N; ++j) { 
        for (k=0; k<mmem_sp3b; k++) mem_sp3b[j][k]=-1; 
        for (k=0; k<mmem_sp3c; k++) mem_sp3c[j][k]=-1; 
        for (k=0; k<mmem_sp4b; k++) mem_sp4b[j][k]=-1; 
        for (k=0; k<mmem_sp4c; k++) mem_sp4c[j][k]=-1;
        for (k=0; k<mmem_sp5b; k++) mem_sp5b[j][k]=-1; 
        for (k=0; k<mmem_sp5c; k++) mem_sp5c[j][k]=-1;
        nmem_sp3b[j]=0;
        nmem_sp3c[j]=0;
        nmem_sp4b[j]=0;
        nmem_sp4c[j]=0;
        nmem_sp5b[j]=0;
        nmem_sp5c[j]=0;
    }

    for (i=0; i<msp4a; ++i) {
        for (j=0;j<4;++j) sp4a[i][j]=-1;
    }
    for (i=0; i<msp4b; ++i) {
        for (j=0;j<5;++j) sp4b[i][j]=-1;
    }
    for (i=0; i<msp4c; ++i) {
        for (j=0;j<6;++j) sp4c[i][j]=-1;
    }
    for (i=0; i<msp5a; ++i) {
        for (j=0;j<5;++j) {
            sp5a[i][j]=-1;
        }
    }
    for (i=0; i<msp5b; ++i) {
        for (j=0;j<6;++j) sp5b[i][j]=-1;
    }
    for (i=0; i<msp5c; ++i) {
        for (j=0;j<7;++j) sp5c[i][j]=-1;
    }
    for (i=0; i<m6Z; ++i) {
        for (j=0;j<6;++j) hc6Z[i][j]=-1;
    }
    
    for (i=0; i<m7K; ++i) {
        for (j=0;j<7;++j) hc7K[i][j]=-1;
    }
    
    for (i=0; i<m8A; ++i) {
        for (j=0;j<8;++j) hc8A[i][j]=-1;
    }
    for (i=0; i<m8B; ++i) {
        for (j=0;j<8;++j) hc8B[i][j]=-1;
    }
    for (i=0; i<m8K; ++i) {
        for (j=0;j<8;++j) hc8K[i][j]=-1;
    }
    
    for (i=0; i<m9A; ++i) {
        for (j=0;j<9;++j) hc9A[i][j]=-1;
    }
    for (i=0; i<m9B; ++i) {
        for (j=0;j<9;++j) hc9B[i][j]=-1;
    }
    for (i=0; i<m9K; ++i) {
        for (j=0;j<9;++j) hc9K[i][j]=-1;
    }
    
    for (i=0; i<m10A; ++i) {
        for (j=0;j<10;++j) hc10A[i][j]=-1;
    }
    for (i=0; i<m10B; ++i) {
        for (j=0;j<10;++j) hc10B[i][j]=-1;
    }
    for (i=0; i<m10K; ++i) {
        for (j=0;j<10;++j) hc10K[i][j]=-1;
    }
    for (i=0; i<m10W; ++i) {
        for (j=0;j<10;++j) hc10W[i][j]=-1;
    }
    
    for (i=0; i<m11A; ++i) {
        for (j=0;j<11;++j) hc11A[i][j]=-1;
    }
    for (i=0; i<m11B; ++i) {
        for (j=0;j<11;++j) hc11B[i][j]=-1;
    }
    for (i=0; i<m11C; ++i) {
        for (j=0;j<11;++j) hc11C[i][j]=-1;
    }
    for (i=0; i<m11E; ++i) {
        for (j=0;j<11;++j) hc11E[i][j]=-1;
    }
    for (i=0; i<m11F; ++i) {
        for (j=0;j<11;++j) hc11F[i][j]=-1;
    }
    for (i=0; i<m11W; ++i) {
        for (j=0;j<11;++j) hc11W[i][j]=-1;
    }
    
    for (i=0; i<m12A; ++i) {
        for (j=0;j<12;++j) hc12A[i][j]=-1;
    }
    for (i=0; i<m12B; ++i) {
        for (j=0;j<12;++j) hc12B[i][j]=-1;
    }
    for (i=0; i<m12D; ++i) {
        for (j=0;j<12;++j) hc12D[i][j]=-1;
    }
    for (i=0; i<m12E; ++i) {
        for (j=0;j<12;++j) hc12E[i][j]=-1;
    }
    for (i=0; i<m12K; ++i) {
        for (j=0;j<12;++j) hc12K[i][j]=-1;
    }
    
    for (i=0; i<m13A; ++i) {
        for (j=0;j<13;++j) hc13A[i][j]=-1;
    }
    for (i=0; i<m13B; ++i) {
        for (j=0;j<13;++j) hc13B[i][j]=-1;
    }
    for (i=0; i<m13K; ++i) {
        for (j=0;j<13;++j) hc13K[i][j]=-1;
    }
    
    for (i=0; i<mFCC; ++i) {
        for (j=0;j<13;++j) hcFCC[i][j]=-1;
    }
    for (i=0; i<mHCP; ++i) {
        for (j=0;j<13;++j) hcHCP[i][j]=-1;
    }
    for (i=0; i<mBCC_9; ++i) {
        for (j=0;j<9;++j) hcBCC_9[i][j]=-1;
    }
    for (i=0; i<mBCC_15; ++i) {
        for (j=0;j<15;++j) hcBCC_15[i][j]=-1;
    }
    
    for (i=0; i<N; ++i) {
        ssp3[i]=ssp3a[i]=ssp3b[i]=ssp3c[i]='C';
        ssp4[i]=ssp4a[i]=ssp4b[i]=ssp4c[i]='C';
        s6Z[i]=s7K[i]='C';
        ssp5[i]=ssp5a[i]=ssp5b[i]=ssp5c[i]='C';
        s8A[i]=s8B[i]=s8K[i]='C';
        s9A[i]=s9B[i]=s9K[i]='C';
        s10A[i]=s10B[i]=s10K[i]=s10W[i]='C';
        s11A[i]=s11B[i]=s11C[i]=s11E[i]=s11F[i]=s11W[i]='C';
        s12A[i]=s12B[i]=s12D[i]=s12E[i]=s12K[i]='C';
        s13A[i]=s13B[i]=s13K[i]='C';
        sFCC[i]=sHCP[i]=sBCC_9[i]=sBCC_15[i]='C';
        
        s9B_cen[i]=s9K_cen[i]='C';
        s10B_cen[i]=s10K_cen[i]=s10W_cen[i]='C';
        s11A_cen[i]=s11B_cen[i]=s11C_cen[i]=s11W_cen[i]='C';
        s12A_cen[i]=s12B_cen[i]=s12K_cen[i]='C';
        s13A_cen[i]=s13B_cen[i]=s13K_cen[i]='C';
        sFCC_cen[i]=sHCP_cen[i]=sBCC_9_cen[i]=sBCC_15_cen[i]='C';
        
        s9B_shell[i]=s9K_shell[i]='C';
        s10B_shell[i]=s10K_shell[i]=s10W_shell[i]='C';
        s11A_shell[i]=s11B_shell[i]=s11C_shell[i]=s11W_shell[i]='C';
        s12A_shell[i]=s12B_shell[i]=s12K_shell[i]='C';
        s13A_shell[i]=s13B_shell[i]=s13K_shell[i]='C';
        sFCC_shell[i]=sHCP_shell[i]=sBCC_9_shell[i]=sBCC_15_shell[i]='C';
    }
}

void Setup_FreeStaticVars()  {  // Free bond detection variables
    int i;
    
    free(rtype);
    free(fXmolName);
    free(fBoxSizeName);
    free(x); free(y); free(z);

    for (i=0; i<N; ++i) {
        free(bNums[i]); 
        free(bondlengths[i]);
    }
    free(bNums); free(bondlengths); free(cnb);
    
    if (USELIST==1) {
        free(map);
        free(head);
        free(llist);
    }

    free(cluster_file_pointers);
    free(raw_file_pointers);


    free(nsp3c_spindlebonds); free(nsp4c_spindlebonds); free(nsp5c_spindlebonds);
    free(nsp3_excess_spindles); free(nsp4_excess_spindles); free(nsp5_excess_spindles);
    
    free(nsp3); free(nsp3a); free(nsp3b); free(nsp3c);
    free(nsp4); free(nsp4a); free(nsp4b); free(nsp4c);
    free(nsp5); free(nsp5a); free(nsp5b); free(nsp5c);
    free(n6Z); free(n7K);
    free(n8A); free(n8B); free(n8K);
    free(n9A); free(n9B); free(n9K);
    free(n10A); free(n10B); free(n10K); free(n10W);
    free(n11A); free(n11B); free(n11C); free(n11E); free(n11F), free(n11W);
    free(n12A); free(n12B); free(n12D); free(n12E); free(n12K);
    free(n13A); free(n13B); free(n13K);
    free(nFCC); free(nHCP); free(nBCC_9); free(nBCC_15);
    
    free(pop_per_frame_sp3); free(pop_per_frame_sp3a); free(pop_per_frame_sp3b); free(pop_per_frame_sp3c);
    free(pop_per_frame_sp4); free(pop_per_frame_sp4a); free(pop_per_frame_sp4b);
    free(pop_per_frame_sp5); free(pop_per_frame_sp5a); free(pop_per_frame_sp5b); free(pop_per_frame_sp5c);
    free(pop_per_frame_sp4c); free(pop_per_frame_6Z); free(pop_per_frame_7K);
    free(pop_per_frame_8A); free(pop_per_frame_8B); free(pop_per_frame_8K);
    free(pop_per_frame_9A); free(pop_per_frame_9B); free(pop_per_frame_9K);
    free(pop_per_frame_10A); free(pop_per_frame_10B); free(pop_per_frame_10K); free(pop_per_frame_10W);
    free(pop_per_frame_11A); free(pop_per_frame_11B); free(pop_per_frame_11C); free(pop_per_frame_11E); free(pop_per_frame_11F); free(pop_per_frame_11W);
    free(pop_per_frame_12A); free(pop_per_frame_12B); free(pop_per_frame_12D); free(pop_per_frame_12E); free(pop_per_frame_12K);
    free(pop_per_frame_13A); free(pop_per_frame_13B); free(pop_per_frame_13K);
    free(pop_per_frame_FCC); free(pop_per_frame_HCP); free(pop_per_frame_BCC_9); free(pop_per_frame_BCC_15);

    for (i=0; i<msp3a; ++i) free(sp3a[i]);
    for (i=0; i<msp3b; ++i) free(sp3b[i]);
    for (i=0; i<msp3c; ++i) free(sp3c[i]);
    for (i=0; i<msp4a; ++i) free(sp4a[i]);
    for (i=0; i<msp4b; ++i) free(sp4b[i]);
    for (i=0; i<msp4c; ++i) free(sp4c[i]);
    for (i=0; i<msp5a; ++i) free(sp5a[i]);
    for (i=0; i<msp5b; ++i) free(sp5b[i]);
    for (i=0; i<msp5c; ++i) free(sp5c[i]);
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
    
    for (i=0; i<N; ++i) free(mem_sp3b[i]);
    for (i=0; i<N; ++i) free(mem_sp3c[i]);
    for (i=0; i<N; ++i) free(mem_sp4b[i]);
    for (i=0; i<N; ++i) free(mem_sp4c[i]);
    for (i=0; i<N; ++i) free(mem_sp5b[i]);
    for (i=0; i<N; ++i) free(mem_sp5c[i]);
    
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

    free(sp3a); free(sp3b); free(sp3c);
    free(sp4a); free(sp4b); free(sp4c);
    free(sp5a); free(sp5b); free(sp5c);
    free(hc6Z); free(hc7K);
    free(hc8A); free(hc8B); free(hc8K);
    free(hc9A); free(hc9B); free(hc9K);
    free(hc10A); free(hc10B); free(hc10K); free(hc10W);
    free(hc11A); free(hc11B); free(hc11C); free(hc11E); free(hc11F); free(hc11W);
    free(hc12A); free(hc12B); free(hc12D); free(hc12E); free(hc12K);
    free(hc13A); free(hc13B); free(hc13K);
    free(hcFCC); free(hcHCP); free(hcBCC_9); free(hcBCC_15);

    free(ssp3); free(ssp3a); free(ssp3b); free(ssp3c);
    free(ssp4); free(ssp4a); free(ssp4b); free(ssp4c);
    free(s6Z); free(s7K);
    free(ssp5); free(ssp5a); free(ssp5b); free(ssp5c);
    free(s8A); free(s8B); free(s8K);
    free(s9A); free(s9B); free(s9K);
    free(s10A); free(s10B); free(s10K); free(s10W);
    free(s11A); free(s11B); free(s11C); free(s11E); free(s11F); free(s11W);
    free(s12A); free(s12B); free(s12D); free(s12E); free(s12K);
    free(s13A); free(s13B); free(s13K);
    free(sFCC);free(sHCP); free(sBCC_9); free(sBCC_15);
    
    free(s9B_cen); free(s9K_cen);
    free(s10B_cen); free(s10K_cen); free(s10W_cen);
    free(s11A_cen); free(s11B_cen); free(s11C_cen); free(s11W_cen); 
    free(s12A_cen); free(s12B_cen); free(s12K_cen);
    free(s13A_cen); free(s13B_cen); free(s13K_cen);
    free(sFCC_cen); free(sHCP_cen); free(sBCC_9_cen); free(sBCC_15_cen);
    
    free(s9B_shell); free(s9K_shell);
    free(s10B_shell); free(s10K_shell); free(s10W_shell);
    free(s11A_shell); free(s11B_shell); free(s11C_shell); free(s11W_shell); 
    free(s12A_shell); free(s12B_shell); free(s12K_shell);
    free(s13A_shell); free(s13B_shell); free(s13K_shell);
    free(sFCC_shell); free(sHCP_shell); free(sBCC_9_shell); free(sBCC_15_shell);
}

void Setup_Centers_Files() {
    char output[1000];

    if (do11AcenXmol==1) {
        printf("\ninitializing 11A centre xmol files...");
        sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_11A_cen.xmol",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
        file_11A_cen_xmol=fopen(output, "w");
        printf("completed\n");
    }

    if (do13AcenXmol==1) {
        printf("\ninitializing 13A centre xmol files...");
        sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_13A_cen.xmol",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
        file_13A_cen_xmol=fopen(output, "w");
        printf("completed\n");
    }
}

void Close_Centers_Files() {
    if (do11AcenXmol==1) {
        printf("Closing 11A centre xmol files....");
        fclose(file_11A_cen_xmol);
        printf("closed!\n\n");
    }

    if (do13AcenXmol==1) {
        printf("Closing 13A centre xmol files....");
        fclose(file_13A_cen_xmol);
        printf("closed!\n\n");
    }
}

int icell(int tix, int tiy, int tiz) { 	// returns cell number (from 1 to ncells) for given (tix,tiy,tiz) coordinate
    return 1 + (tix-1+M)%M + M*((tiy-1+M)%M) + M*M*((tiz-1+M)%M);
}

void Setup_Cell_List() {

    int ix, iy, iz;
    int imap;

    if (USELIST==1) {
        M = (int)(side/rcutAA);	// number of cells along box side
        if (M<3) Error_no_free("main(): M<3, too few cells");
        ncells = M*M*M;	// total number of cells
        cellSide = side/M;	// length of cells
        invcellSide = 1.0/cellSide;	// invcellSide
        printf("m %d ncells %d cellside %.15lg\n", M, ncells, cellSide);
        // routine to create the thirteen nearest neighbours array map[] of each cell
        for (iz=1; iz<=M; iz++) {
            for (iy=1; iy<=M; iy++) {
                for (ix=1; ix<=M; ix++) {
                    imap = (icell(ix,iy,iz)-1)*13;
                    map[imap+1 ]=icell(ix+1,iy	,iz	);
                    map[imap+2 ]=icell(ix+1,iy+1,iz	);
                    map[imap+3 ]=icell(ix	 ,iy+1,iz	);
                    map[imap+4 ]=icell(ix-1 ,iy+1,iz	);
                    map[imap+5 ]=icell(ix+1,iy	,iz-1	);
                    map[imap+6 ]=icell(ix+1,iy+1,iz-1	);
                    map[imap+7 ]=icell(ix	 ,iy+1,iz-1	);
                    map[imap+8 ]=icell(ix-1 ,iy+1,iz-1	);
                    map[imap+9 ]=icell(ix+1,iy	,iz+1	);
                    map[imap+10]=icell(ix+1,iy+1,iz+1	);
                    map[imap+11]=icell(ix	 ,iy+1,iz+1	);
                    map[imap+12]=icell(ix-1 ,iy+1,iz+1);
                    map[imap+13]=icell(ix	 ,iy	,iz+1	);
                }
            }
        }
    }
}