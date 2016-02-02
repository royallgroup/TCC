#include "setup.h"
#include "tools.h"
#include "iniparser.h"

//// START: Setup routines
int Setup_GetFirstIntFromLine(FILE *stream) {   // take and return first integer from line in file stream
    char input[10000], errMsg[1000];
    char *pch;
    if (fgets(input,10000,stream)==NULL) {
        sprintf(errMsg,"Setup_GetFirstIntFromLine(): end of input file reached\n");
        Error_no_free(errMsg);
    }
    pch=strtok(input," ");

    return atoi(pch);
}

double Setup_GetFirstDoubleFromLine(FILE *stream) { // take and return first double from line in file stream
    char input[10000], errMsg[1000];
    char *pch;
    if (fgets(input,10000,stream)==NULL) {
        sprintf(errMsg,"Setup_GetFirstDoubleFromLine(): end of input file reached\n");
        Error_no_free(errMsg);
    }
    pch=strtok(input," ");
    return atof(pch);
}

void Setup_ReadIniFile(char *filename) {
    
    char errMsg[1000];
    double RHO;
    dictionary  *   ini ;

    fXmolName=malloc(500*sizeof(char)); if (fXmolName==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): fXmolName[] malloc out of memory\n");   Error_no_free(errMsg); }
    fBoxSizeName=malloc(500*sizeof(char)); if (fXmolName==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): fBoxSizeName[] malloc out of memory\n");   Error_no_free(errMsg); }

    ini = iniparser_load(filename);
    if (ini==NULL) {
        sprintf(errMsg,"Setup_ReadIniFile(): Error opening file %s",filename);  // Always test file open
        Error_no_free(errMsg);
    }
    
    //box
    ISNOTCUBIC = iniparser_getint(ini, "box:box_type", -1);
    fBoxSizeName = (char*)iniparser_getstring(ini, "box:box_name", "-1");
    
    
    //run
    fXmolName = (char*)iniparser_getstring(ini, "run:xyzfilename", "-1");
    FRAMES = iniparser_getint(ini, "run:frames", -1);
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
    WHICHPOTENTIAL = iniparser_getint(ini, "simulation:potential_type", -1);

    //output
    doWriteBonds = iniparser_getboolean(ini, "output:bonds", -1);
    doWriteClus = iniparser_getboolean(ini, "output:clusts", -1);
    doWriteRaw = iniparser_getboolean(ini, "output:raw", -1);
    do11AcenXmol = iniparser_getboolean(ini, "output:11a", -1);
    do13AcenXmol = iniparser_getboolean(ini, "output:13a", -1);
    doWritePopPerFrame = iniparser_getboolean(ini, "output:pop_per_frame", -1);
    binWidth = iniparser_getdouble(ini, "output:bin_width", -1);
    doBLDistros = iniparser_getboolean(ini, "output:bond_length", -1);
    doClusBLDistros = iniparser_getboolean(ini, "output:bond_length_cluster", -1);
    doClusBLDeviation = iniparser_getboolean(ini, "output:bond_length_dev", -1);
    donbDistros = iniparser_getboolean(ini, "output:neighbour_dist", -1);
    doBondedCen = iniparser_getboolean(ini, "output:bonded_dist", -1);
    doClusComp = iniparser_getboolean(ini, "output:cluster_composition", -1);
    doSubClusts = iniparser_getboolean(ini, "output:subclusters", -1);

    doPotential = iniparser_getboolean(ini, "extra:potential_energy", -1);
    doCoslovich = iniparser_getboolean(ini, "extra:coslovich", -1);
    doDynamics = iniparser_getboolean(ini, "extra:dodynamics", -1);
    talpha = iniparser_getdouble(ini, "extra:alpha_time", -1);
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
    initNoLifetimes=initNoDynamicClusters=incrDynamicClusters=1;
    if (doDynamics==1) {
        printf("In dynamic mode therefore not writing any files other than those for the dynamics of clusters\n");
        doWriteBonds=doWriteClus=doWriteRaw=doWritePopPerFrame=doBLDistros=doClusBLDistros=doClusBLDeviation=donbDistros=doBondedCen=doClusComp=doPotential=doCoslovich=0;
    }
    
    if (NA>N) {
        Error_no_free("Setup_ReadIniFile(): NA > N - something wrong in xmol trajectory params file\n");
    }
    else if (NA<N) {doBinary=1;}
    else {doBinary=0;}
    
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
    printf("binWidth %.5lg doBLDistros %d doClusBLDistros %d doClusBLDeviation %d\n",binWidth,doBLDistros,doClusBLDistros,doClusBLDeviation);
    printf("donbDistros %d doBondedCen %d doClusComp %d doSubClusts %d\n",donbDistros,doBondedCen,doClusComp,doSubClusts);
    printf("doPotential %d doCoslovich %d doDynamics %d talpha %lg PRINTINFO %d\n\n",doPotential,doCoslovich,doDynamics,talpha,PRINTINFO);
        
    if (ISNOTCUBIC==0) {
        printf("calculating box sides from RHO\n");
        printf("box side length = %.5lg, half side: %.5lg\n\n",sidex,halfSidex);
    }
    
    if (doPotential==1) {
        if (WHICHPOTENTIAL==0) printf("%d: (Binary) Lennard-Jones potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==1) printf("%d: Stoddard-Ford cut-off (Binary) Lennard-Jones potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==2) printf("%d: Morse+Yukawa potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==3) {
            printf("%d: Not currently implemented potential\n",WHICHPOTENTIAL);
            sprintf(errMsg,"Setup_ReadIniFile(): Not currently implemented potential"); // Always test file open
            Error_no_free(errMsg);
        }
        else if (WHICHPOTENTIAL==4) printf("%d: Inverse-power law potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==5) printf("%d: (Binary) WCA_s cubic splined potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==6) printf("%d: Stoddard-Ford cut-off (Binary) Inverse-power law potential\n",WHICHPOTENTIAL);
        else if (WHICHPOTENTIAL==7) printf("%d: CRVT potential\n",WHICHPOTENTIAL);
        else {
            sprintf(errMsg,"Setup_ReadIniFile(): WHICHPOTENTIAL %d not implemented",WHICHPOTENTIAL);    // Always test file open
            Error_no_free(errMsg);
        }
    }
    iniparser_freedict(ini);
}

void Setup_ReadBox(FILE *readIn)  {
    int sweep;
        
    if (feof(readIn)) Error("Setup_ReadBox(): end of input file reached\n");
    if(ISNOTCUBIC==3){
        printf("Triclinic Boundary Conditions\n");
      
        fscanf(readIn,"%i %lg %lg %lg %lg\n", &sweep,&sidex,&sidey,&sidez, &tilt);
    
        printf("%i %lg %lg %lg %lg\n",sweep,sidex,sidey,sidez, tilt);}
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
        if (PRINTINFO==1) printf("d%d Setup_Readxyz(): TCC analysis frame %d/%d - reading in %d particles frame %d of XMOL\n",rank,f,FRAMES-1,i,e);
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
            if (PRINTINFO==1) if (i==N-1) printf("d%d f%d part%d %c %.5lg %.5lg %.5lg\n\n",rank,f,i,c,x[i],y[i],z[i]);
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
    double lengthstemp;

    dosp3=dosp3a=dosp3b=dosp3c=1;
    dosp4=dosp4a=dosp4b=dosp4c=1;
    dosp5=dosp5a=dosp5b=dosp5c=1;
    do6Z=do7K=do8A=do8B=do8K=do9A=do9B=do9K=do10A=do10B=do10K=do10W=1;
    do11A=do11B=do11C=do11E=do11F=do11W=do12A=do12B=do12D=do12E=do12K=1;
    do13A=do13B=do13K=doFCC=doHCP=doBCC9=doBCC15=1;
    
    msp3a=msp3b=msp3c=initNoStatic; // max size of **sp** arrays in dimension i
    msp4a=msp4b=msp4c=m6A=initNoStatic; // max size of **sp** arrays in dimension i
    msp5a=msp5b=msp5c=initNoStatic; // max size of **sp** arrays in dimension i
    m6Z=m7K=initNoStatic;   // max size of m** arrays in dimension i
    m8A=m8B=m8K=initNoStatic;   // max size of m** arrays in dimension i
    m9A=m9B=m9K=initNoStatic;   // max size of m** arrays in dimension i
    m10A=m10B=m10K=m10W=initNoStatic;   // max size of m** arrays in dimension i
    m11A=m11B=m11C=m11E=m11F=m11W=initNoStatic; // max size of m** arrays in dimension i
    m12A=m12B=m12D=m12E=m12K=initNoStatic;  // max size of m** arrays in dimension i
    m13A=m13B=m13K=initNoStatic;    // max size of m** arrays in dimension i
    mFCC=mHCP=mBCC_9=mBCC_15=initNoStatic;  // max size of **sp** arrays in dimension i
    
    msp3=msp3a+msp3b+msp3c;
    msp4=msp4a+msp4b+msp4c;
    msp5=msp5a+msp5b+msp5c;
    
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
    mean_pop_per_frame_sp4=mean_pop_per_frame_sp4a=mean_pop_per_frame_sp4b=mean_pop_per_frame_6A=0.0;
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

    maxnb=totNclus=0;
    maxto3=maxto4=maxto5=0;
    
    // count number of bonds between spindle particles
    nsp3c_spindlebonds = malloc(FRAMES*sizeof(int));    if (nsp3c_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp3c_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
    nsp4c_spindlebonds = malloc(FRAMES*sizeof(int));    if (nsp4c_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nsp4c_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
    n6A_spindlebonds = malloc(FRAMES*sizeof(int));  if (n6A_spindlebonds==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n6A_spindlebonds[] malloc out of memory\n");  Error_no_free(errMsg); }
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
    n6A = malloc(FRAMES*sizeof(int));   if (n6A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n6A[] malloc out of memory\n");    Error_no_free(errMsg); }
    
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
        n6A_spindlebonds[f]=0;
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
        n6A[f]=0;
        
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
    hc6A = malloc(m6A*sizeof(int *));   if (hc6A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc6A[] malloc out of memory\n");  Error_no_free(errMsg); }
    for (j=0; j<msp4c; ++j) { hc6A[j] = malloc(6*sizeof(int));  if (hc6A[j]==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): hc6A[][] malloc out of memory\n"); Error_no_free(errMsg); } }
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
    for (j=0; j<msp4c; ++j) { for (k=0; k<6; k++) { sp4c[j][k]=-1; hc6A[j][k]=-1; }}
        
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
    s5A=malloc(N*sizeof(char)); if (s5A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s5A[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp4=malloc(N*sizeof(char)); if (ssp4==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4a=malloc(N*sizeof(char)); if (ssp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp4b=malloc(N*sizeof(char)); if (ssp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp4b[] malloc out of memory\n"); Error_no_free(errMsg); }
    s6A=malloc(N*sizeof(char)); if (s6A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s6A[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s6Z=malloc(N*sizeof(char)); if (s6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s6Z[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    s7K=malloc(N*sizeof(char)); if (s7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s7K[] malloc out of memory\n"); Error_no_free(errMsg); }
    
    ssp5=malloc(N*sizeof(char)); if (ssp5==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5a=malloc(N*sizeof(char)); if (ssp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5a[] malloc out of memory\n"); Error_no_free(errMsg); }
    ssp5b=malloc(N*sizeof(char)); if (ssp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): ssp5b[] malloc out of memory\n"); Error_no_free(errMsg); }
    s7A=malloc(N*sizeof(char)); if (s7A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): s7A[] malloc out of memory\n"); Error_no_free(errMsg); }
    
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
        s5A[j]='C';
        
        ssp4[j]='C';
        ssp4a[j]='C';
        ssp4b[j]='C';
        s6A[j]='C';
        
        ssp5[j]='C';
        ssp5a[j]='C';
        ssp5b[j]='C';
        s7A[j]='C';
        
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
    pop_per_frame_6A = malloc(FRAMES*sizeof(double));   if (pop_per_frame_6A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): pop_per_frame_6A[] malloc out of memory\n");  Error_no_free(errMsg); }
    
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
        pop_per_frame_6A[f]=0;
        
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
    
    if (doBLDistros) { // histograms of bond lengths
        lengthstemp=rcutAA; // calculate number of bins for histogram of bond lengths
        if (rcutAB>lengthstemp) lengthstemp=rcutAB;
        if (rcutBB>lengthstemp) lengthstemp=rcutBB;
        BLDistroNoBins=(int)ceil(lengthstemp/binWidth);
        
        BLDistroNoSamples=0;
        meanBL=0.0;
        BLDistro=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro[j]=0;
        if (doBinary==1) {
            BLDistroNoSamplesAA=BLDistroNoSamplesAB=BLDistroNoSamplesBB=0;
            meanBLAA=meanBLAB=meanBLBB=0.0;
            
            BLDistroAA=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroAA==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroAA[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<BLDistroNoBins; j++) BLDistroAA[j]=0;
        
            BLDistroAB=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroAB==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroAB[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<BLDistroNoBins; j++) BLDistroAB[j]=0;

            BLDistroBB=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroBB==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroBB[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<BLDistroNoBins; j++) BLDistroBB[j]=0;
        }
    }
    
    if (donbDistros==1) { // distributions of number of bonds to each particle
        nbDistroNoSamples=0;
        meannb=0.0;
        nbDistro=malloc((nB+1)*sizeof(int)); if (nbDistro==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nbDistro[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<(nB+1); j++) nbDistro[j]=0;
        
        if (doBinary==1) {
            nbDistroNoSamplesAA=nbDistroNoSamplesAB=nbDistroNoSamplesBA=nbDistroNoSamplesBB=0;  // look at statistics of bonds (mean number of bonds per part, mean number AA bonds per A part etc etc)
            meannbAA=meannbAB=meannbBA=meannbBB=0.0;
            
            nbDistroAA=malloc((nB+1)*sizeof(int)); if (nbDistroAA==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nbDistroAA[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<(nB+1); j++) nbDistroAA[j]=0;
        
            nbDistroAB=malloc((nB+1)*sizeof(int)); if (nbDistroAB==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nbDistroAB[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<(nB+1); j++) nbDistroAB[j]=0;
        
            nbDistroBA=malloc((nB+1)*sizeof(int)); if (nbDistroBA==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nbDistroBA[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<(nB+1); j++) nbDistroBA[j]=0;
        
            nbDistroBB=malloc((nB+1)*sizeof(int)); if (nbDistroBB==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): nbDistroBB[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<(nB+1); j++) nbDistroBB[j]=0;
        }
    }
    
    if (doBondedCen==1) {   // distribution of number of particles bonded to central particle of a cluster
        n_distro_bonded_to_cen_9B=malloc((nB+1)*sizeof(int));   if (n_distro_bonded_to_cen_9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_9B[] malloc out of memory\n");    Error_no_free(errMsg); }
        n_distro_bonded_to_cen_9K=malloc((nB+1)*sizeof(int));   if (n_distro_bonded_to_cen_9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_9K[] malloc out of memory\n");    Error_no_free(errMsg); }
        n_distro_bonded_to_cen_10B=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_10B[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_10K=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_10K[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_10W=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_10W[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_11A=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_11A[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_11B=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_11B[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_11C=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_11C[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_11W=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_11W[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_12A=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_12A[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_12B=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_12B[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_12K=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_12K[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_13A=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_13A[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_13K=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_13K[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_13B=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_13B[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_FCC=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_FCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_FCC[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_HCP=malloc((nB+1)*sizeof(int));  if (n_distro_bonded_to_cen_HCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_HCP[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_BCC_9=malloc((nB+1)*sizeof(int));    if (n_distro_bonded_to_cen_BCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_BCC_9[] malloc out of memory\n");  Error_no_free(errMsg); }
        n_distro_bonded_to_cen_BCC_15=malloc((nB+1)*sizeof(int));   if (n_distro_bonded_to_cen_BCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): n_distro_bonded_to_cen_BCC_15[] malloc out of memory\n");    Error_no_free(errMsg); }
        
        for (j=0; j<nB+1; j++) {    // set variables to zero
            n_distro_bonded_to_cen_9B[j]=0;
            n_distro_bonded_to_cen_9K[j]=0;
            n_distro_bonded_to_cen_10B[j]=0;
            n_distro_bonded_to_cen_10K[j]=0;
            n_distro_bonded_to_cen_10W[j]=0;
            n_distro_bonded_to_cen_11A[j]=0;
            n_distro_bonded_to_cen_11B[j]=0;
            n_distro_bonded_to_cen_11C[j]=0;
            n_distro_bonded_to_cen_11W[j]=0;
            n_distro_bonded_to_cen_12A[j]=0;
            n_distro_bonded_to_cen_12B[j]=0;
            n_distro_bonded_to_cen_12K[j]=0;
            n_distro_bonded_to_cen_13A[j]=0;
            n_distro_bonded_to_cen_13B[j]=0;
            n_distro_bonded_to_cen_13K[j]=0;
            n_distro_bonded_to_cen_FCC[j]=0;
            n_distro_bonded_to_cen_HCP[j]=0;
            n_distro_bonded_to_cen_BCC_9[j]=0;
            n_distro_bonded_to_cen_BCC_15[j]=0;
        }
        // number of samples for number of particles bonded to central particle distributions
        n_bonded_to_cen_9B=n_bonded_to_cen_9K=0;
        n_bonded_to_cen_10B=n_bonded_to_cen_10K=n_bonded_to_cen_10W=0;
        n_bonded_to_cen_11A=n_bonded_to_cen_11B=n_bonded_to_cen_11C=n_bonded_to_cen_11W=0;
        n_bonded_to_cen_12A=n_bonded_to_cen_12B=n_bonded_to_cen_12K=0;
        n_bonded_to_cen_13A=n_bonded_to_cen_13B=n_bonded_to_cen_13K=0;
        n_bonded_to_cen_FCC=n_bonded_to_cen_HCP=n_bonded_to_cen_BCC_9=n_bonded_to_cen_BCC_15=0;
    }
    
    if (doClusBLDistros==1) {   // do bond length distributions between neighbouring particles within a particular cluster
        BLDistroNoSamplessp3=BLDistroNoSamplessp3a=BLDistroNoSamplessp3b=BLDistroNoSamplessp3c=0;   // number of samples
        BLDistroNoSamplessp4=BLDistroNoSamplessp4a=BLDistroNoSamplessp4b=BLDistroNoSamplessp4c=0;
        BLDistroNoSamples6A=0;
        BLDistroNoSamplessp5=BLDistroNoSamplessp5a=BLDistroNoSamplessp5b=BLDistroNoSamplessp5c=0;
        BLDistroNoSamples6Z=BLDistroNoSamples7K=0;
        BLDistroNoSamples8A=BLDistroNoSamples8B=BLDistroNoSamples8K=0;  
        BLDistroNoSamples9A=BLDistroNoSamples9B=BLDistroNoSamples9K=0;
        BLDistroNoSamples10A=BLDistroNoSamples10B=BLDistroNoSamples10K=BLDistroNoSamples10W=0;
        BLDistroNoSamples11A=BLDistroNoSamples11B=BLDistroNoSamples11C=BLDistroNoSamples11E=BLDistroNoSamples11F=BLDistroNoSamples11W=0;
        BLDistroNoSamples12A=BLDistroNoSamples12B=BLDistroNoSamples12D=BLDistroNoSamples12E=BLDistroNoSamples12K=0;
        BLDistroNoSamples13A=BLDistroNoSamples13B=BLDistroNoSamples13K=0;
        BLDistroNoSamplesFCC=BLDistroNoSamplesHCP=BLDistroNoSamplesBCC_9=BLDistroNoSamplesBCC_15=0;
        
        meanBLsp3=meanBLsp3a=meanBLsp3b=meanBLsp3c=0;   // mean bond length
        meanBLsp4=meanBLsp4a=meanBLsp4b=meanBLsp4c=0;
        meanBL6A=0;
        meanBLsp5=meanBLsp5a=meanBLsp5b=meanBLsp5c=0;
        meanBL6Z=meanBL7K=0;
        meanBL8A=meanBL8B=meanBL8K=0;   
        meanBL9A=meanBL9B=meanBL9K=0;
        meanBL10A=meanBL10B=meanBL10K=meanBL10W=0;
        meanBL11A=meanBL11B=meanBL11C=meanBL11E=meanBL11F=meanBL11W=0;
        meanBL12A=meanBL12B=meanBL12D=meanBL12E=meanBL12K=0;
        meanBL13A=meanBL13B=meanBL13K=0;
        meanBLFCC=meanBLHCP=meanBLBCC_9=meanBLBCC_15=0;

        BLDistrosp3=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp3==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp3[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp3[j]=0;
        BLDistrosp3a=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp3a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp3a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp3a[j]=0;
        BLDistrosp3b=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp3b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp3b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp3b[j]=0;
        BLDistrosp3c=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp3c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp3c[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp3c[j]=0;

        BLDistrosp4=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp4==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp4[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp4[j]=0;
        BLDistrosp4a=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp4a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp4a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp4a[j]=0;
        BLDistrosp4b=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp4b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp4b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp4b[j]=0;
        BLDistrosp4c=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp4c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp4c[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp4c[j]=0;
        BLDistro6A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro6A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro6A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro6A[j]=0;
        
        BLDistrosp5=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp5==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp5[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp5[j]=0;
        BLDistrosp5a=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp5a==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp5a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp5a[j]=0;
        BLDistrosp5b=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp5b==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp5b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp5b[j]=0;
        BLDistrosp5c=malloc(BLDistroNoBins*sizeof(int)); if (BLDistrosp5c==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistrosp5c[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistrosp5c[j]=0;
        BLDistro6Z=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro6Z==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro6Z[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro6Z[j]=0;
        BLDistro7K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro7K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro7K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro7K[j]=0;
        BLDistro8A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro8A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro8A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro8A[j]=0;
        BLDistro8B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro8B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro8B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro8B[j]=0;
        BLDistro8K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro8K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro8K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro8K[j]=0;
        BLDistro9A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro9A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro9A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro9A[j]=0;
        BLDistro9B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro9B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro9B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro9B[j]=0;
        BLDistro9K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro9K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro9K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro9K[j]=0;
        BLDistro10A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro10A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro10A[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro10A[j]=0;
        BLDistro10B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro10B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro10B[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro10B[j]=0;
        BLDistro10K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro10K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro10K[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro10K[j]=0;
        BLDistro10W=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro10W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro10W[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro10W[j]=0;
        BLDistro11A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11A[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11A[j]=0;
        BLDistro11B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11B[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11B[j]=0;
        BLDistro11C=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11C==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11C[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11C[j]=0;
        BLDistro11E=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11E[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11E[j]=0;
        BLDistro11F=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11F==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11F[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11F[j]=0;
        BLDistro11W=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro11W==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro11W[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro11W[j]=0;
        BLDistro12A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro12A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro12A[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro12A[j]=0;
        BLDistro12B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro12B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro12B[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro12B[j]=0;
        BLDistro12D=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro12D==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro12D[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro12D[j]=0;
        BLDistro12E=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro12E==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro12E[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro12E[j]=0;
        BLDistro12K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro12K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro12K[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro12K[j]=0;
        BLDistro13A=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro13A==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro13A[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro13A[j]=0;
        BLDistro13B=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro13B==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro13B[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro13B[j]=0;
        BLDistro13K=malloc(BLDistroNoBins*sizeof(int)); if (BLDistro13K==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistro13K[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistro13K[j]=0;
        BLDistroFCC=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroFCC==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroFCC[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistroFCC[j]=0;
        BLDistroHCP=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroHCP==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroHCP[] malloc out of memory\n");    Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistroHCP[j]=0;
        BLDistroBCC_9=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroBCC_9==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroBCC_9[] malloc out of memory\n");  Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistroBCC_9[j]=0;
        BLDistroBCC_15=malloc(BLDistroNoBins*sizeof(int)); if (BLDistroBCC_15==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): BLDistroBCC_15[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<BLDistroNoBins; j++) BLDistroBCC_15[j]=0;
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
    
    nsp3c_spindlebonds[f]=nsp4c_spindlebonds[f]=n6A_spindlebonds[f]=nsp5c_spindlebonds[f]=0;
    nsp3_excess_spindles[f]=nsp4_excess_spindles[f]=nsp5_excess_spindles[f]=0;
    
    nsp3[f]=nsp3a[f]=nsp3b[f]=nsp3c[f]=0;
    nsp4[f]=nsp4a[f]=nsp4b[f]=nsp4c[f]=0;
    nsp5[f]=nsp5a[f]=nsp5b[f]=nsp5c[f]=0;
    
    n6A[f]=n6Z[f]=n7K[f]=0;
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
    pop_per_frame_6A[f]=pop_per_frame_6Z[f]=pop_per_frame_7K[f]=0.0;
    pop_per_frame_8A[f]=pop_per_frame_8B[f]=pop_per_frame_8K[f]=0.0;
    pop_per_frame_9A[f]=pop_per_frame_9B[f]=pop_per_frame_9K[f]=0.0;
    pop_per_frame_10A[f]=pop_per_frame_10B[f]=pop_per_frame_10K[f]=pop_per_frame_10W[f]=0.0;
    pop_per_frame_11A[f]=pop_per_frame_11B[f]=pop_per_frame_11C[f]=pop_per_frame_11E[f]=pop_per_frame_11F[f]=pop_per_frame_11W[f]=0.0;
    pop_per_frame_12A[f]=pop_per_frame_12B[f]=pop_per_frame_12D[f]=pop_per_frame_12E[f]=pop_per_frame_12K[f]=0.0;
    pop_per_frame_13A[f]=pop_per_frame_13B[f]=pop_per_frame_13K[f]=0.0;
    pop_per_frame_FCC[f]=pop_per_frame_HCP[f]=pop_per_frame_BCC_9[f]=pop_per_frame_BCC_15[f]=0.0;
    
    if (doClusBLDeviation==1) {
        for (i=0; i<msp3; ++i) {
            bl_mom_sp3[i]=0.0;
        }
    }
    for (i=0; i<msp3a; ++i) {
        for (j=0;j<3;++j) sp3a[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp3a[i]=0.0;
    }
    for (i=0; i<msp3b; ++i) {
        for (j=0;j<4;++j) sp3b[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp3b[i]=0.0;
    }
    for (i=0; i<msp3c; ++i) {
        for (j=0;j<5;++j) sp3c[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp3c[i]=0.0;
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
    if (doClusBLDeviation==1) {
        for (i=0; i<msp4; ++i) {
            bl_mom_sp4[i]=0.0;
        }
    }
    for (i=0; i<msp4a; ++i) {
        for (j=0;j<4;++j) sp4a[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp4a[i]=0.0;
    }
    for (i=0; i<msp4b; ++i) {
        for (j=0;j<5;++j) sp4b[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp4b[i]=0.0;
    }
    for (i=0; i<msp4c; ++i) {
        for (j=0;j<6;++j) sp4c[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp4c[i]=0.0;
    }
    for (i=0; i<m6A; ++i) {
        for (j=0;j<6;++j) hc6A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_6A[i]=0.0;
    }
    
    if (doClusBLDeviation==1) {
        for (i=0; i<msp5; ++i) {
            bl_mom_sp5[i]=0.0;
        }
    }
    for (i=0; i<msp5a; ++i) { 
        for (j=0;j<5;++j) {
            sp5a[i][j]=-1;
        }
        if (doClusBLDeviation==1) bl_mom_sp5a[i]=0.0;
    }
    for (i=0; i<msp5b; ++i) {
        for (j=0;j<6;++j) sp5b[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp5b[i]=0.0;
    }
    for (i=0; i<msp5c; ++i) {
        for (j=0;j<7;++j) sp5c[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_sp5c[i]=0.0;
    }
    
    for (i=0; i<m6Z; ++i) {
        for (j=0;j<6;++j) hc6Z[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_6Z[i]=0.0;
    }
    
    for (i=0; i<m7K; ++i) {
        for (j=0;j<7;++j) hc7K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_7K[i]=0.0;
    }
    
    for (i=0; i<m8A; ++i) {
        for (j=0;j<8;++j) hc8A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_8A[i]=0.0;
    }
    for (i=0; i<m8B; ++i) {
        for (j=0;j<8;++j) hc8B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_8B[i]=0.0;
    }
    for (i=0; i<m8K; ++i) {
        for (j=0;j<8;++j) hc8K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_8K[i]=0.0;
    }
    
    for (i=0; i<m9A; ++i) {
        for (j=0;j<9;++j) hc9A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_9A[i]=0.0;
    }
    for (i=0; i<m9B; ++i) {
        for (j=0;j<9;++j) hc9B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_9B[i]=0.0;
    }
    for (i=0; i<m9K; ++i) {
        for (j=0;j<9;++j) hc9K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_9K[i]=0.0;
    }
    
    for (i=0; i<m10A; ++i) {
        for (j=0;j<10;++j) hc10A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_10A[i]=0.0;
    }
    for (i=0; i<m10B; ++i) {
        for (j=0;j<10;++j) hc10B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_10B[i]=0.0;
    }
    for (i=0; i<m10K; ++i) {
        for (j=0;j<10;++j) hc10K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_10K[i]=0.0;
    }
    for (i=0; i<m10W; ++i) {
        for (j=0;j<10;++j) hc10W[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_10W[i]=0.0;
    }
    
    for (i=0; i<m11A; ++i) {
        for (j=0;j<11;++j) hc11A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11A[i]=0.0;
    }
    for (i=0; i<m11B; ++i) {
        for (j=0;j<11;++j) hc11B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11B[i]=0.0;
    }
    for (i=0; i<m11C; ++i) {
        for (j=0;j<11;++j) hc11C[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11C[i]=0.0;
    }
    for (i=0; i<m11E; ++i) {
        for (j=0;j<11;++j) hc11E[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11E[i]=0.0;
    }
    for (i=0; i<m11F; ++i) {
        for (j=0;j<11;++j) hc11F[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11F[i]=0.0;
    }
    for (i=0; i<m11W; ++i) {
        for (j=0;j<11;++j) hc11W[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_11W[i]=0.0;
    }
    
    for (i=0; i<m12A; ++i) {
        for (j=0;j<12;++j) hc12A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_12A[i]=0.0;
    }
    for (i=0; i<m12B; ++i) {
        for (j=0;j<12;++j) hc12B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_12B[i]=0.0;
    }
    for (i=0; i<m12D; ++i) {
        for (j=0;j<12;++j) hc12D[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_12D[i]=0.0;
    }
    for (i=0; i<m12E; ++i) {
        for (j=0;j<12;++j) hc12E[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_12E[i]=0.0;
    }
    for (i=0; i<m12K; ++i) {
        for (j=0;j<12;++j) hc12K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_12K[i]=0.0;
    }
    
    for (i=0; i<m13A; ++i) {
        for (j=0;j<13;++j) hc13A[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_13A[i]=0.0;
    }
    for (i=0; i<m13B; ++i) {
        for (j=0;j<13;++j) hc13B[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_13B[i]=0.0;
    }
    for (i=0; i<m13K; ++i) {
        for (j=0;j<13;++j) hc13K[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_13K[i]=0.0;
    }
    
    for (i=0; i<mFCC; ++i) {
        for (j=0;j<13;++j) hcFCC[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_FCC[i]=0.0;
    }
    for (i=0; i<mHCP; ++i) {
        for (j=0;j<13;++j) hcHCP[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_HCP[i]=0.0;
    }
    for (i=0; i<mBCC_9; ++i) {
        for (j=0;j<9;++j) hcBCC_9[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_BCC_9[i]=0.0;
    }
    for (i=0; i<mBCC_15; ++i) {
        for (j=0;j<15;++j) hcBCC_15[i][j]=-1;
        if (doClusBLDeviation==1) bl_mom_BCC_15[i]=0.0;
    }
    
    for (i=0; i<N; ++i) {
        ssp3[i]=ssp3a[i]=ssp3b[i]=s5A[i]='C';
        ssp4[i]=ssp4a[i]=ssp4b[i]=s6A[i]='C';
        s6Z[i]=s7K[i]='C';
        ssp5[i]=ssp5a[i]=ssp5b[i]=s7A[i]='C';
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
    
    if (doPotential==1) {
        potential=0.0;
        for (i=0; i<N; ++i) part_pot[i]=0.0;
    }
    
    if (doSubClusts==1 && doDynamics==1) {
        if (dyn_msp3b!=-1) for (i=0; i<msp3b; ++i) dyn_up_sp3b[i]=-1;
        if (dyn_msp3c!=-1) for (i=0; i<msp3c; ++i) dyn_up_sp3c[i]=-1;
        if (dyn_msp4b!=-1) for (i=0; i<msp4b; ++i) dyn_up_sp4b[i]=-1;
        if (dyn_m6A!=-1) for (i=0; i<msp4c; ++i) dyn_up_sp4c[i]=-1;
        if (dyn_msp5b!=-1) for (i=0; i<msp5b; ++i) dyn_up_sp5b[i]=-1;
        if (dyn_msp5c!=-1) for (i=0; i<msp5c; ++i) dyn_up_sp5c[i]=-1;
        if (dyn_m9B!=-1) for (i=0; i<m9B; ++i) dyn_up_9B[i]=-1;
        if (dyn_m9K!=-1) for (i=0; i<m9K; ++i) dyn_up_9K[i]=-1;
        if (dyn_m10B!=-1) for (i=0; i<m10B; ++i) dyn_up_10B[i]=-1;
        if (dyn_m11A!=-1) for (i=0; i<m11A; ++i) dyn_up_11A[i]=-1;
        if (dyn_m11C!=-1) for (i=0; i<m11C; ++i) dyn_up_11C[i]=-1;
        if (dyn_m11F!=-1) for (i=0; i<m11F; ++i) dyn_up_11F[i]=-1;
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

    free(nsp3c_spindlebonds); free(nsp4c_spindlebonds); free(n6A_spindlebonds); free(nsp5c_spindlebonds);
    free(nsp3_excess_spindles); free(nsp4_excess_spindles); free(nsp5_excess_spindles);
    
    free(nsp3); free(nsp3a); free(nsp3b); free(nsp3c);
    free(nsp4); free(nsp4a); free(nsp4b); free(nsp4c);
    free(nsp5); free(nsp5a); free(nsp5b); free(nsp5c);
    free(n6A); free(n6Z); free(n7K);
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
    free(pop_per_frame_6A); free(pop_per_frame_6Z); free(pop_per_frame_7K);
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
    for (i=0; i<m6A; ++i) free(hc6A[i]);
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
    free(hc6A); free(hc6Z); free(hc7K);
    free(hc8A); free(hc8B); free(hc8K);
    free(hc9A); free(hc9B); free(hc9K);
    free(hc10A); free(hc10B); free(hc10K); free(hc10W);
    free(hc11A); free(hc11B); free(hc11C); free(hc11E); free(hc11F); free(hc11W);
    free(hc12A); free(hc12B); free(hc12D); free(hc12E); free(hc12K);
    free(hc13A); free(hc13B); free(hc13K);
    free(hcFCC); free(hcHCP); free(hcBCC_9); free(hcBCC_15);

    free(ssp3); free(ssp3a); free(ssp3b); free(s5A);
    free(ssp4); free(ssp4a); free(ssp4b); free(s6A); 
    free(s6Z); free(s7K);
    free(ssp5); free(ssp5a); free(ssp5b); free(s7A);
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

    if (doBLDistros==1) {
        free(BLDistro);
        if (doBinary==1) {
            free(BLDistroAA); free(BLDistroAB); free(BLDistroBB); 
            
        }
    }

    if (donbDistros==1) {
        free(nbDistro); 
        if (doBinary==1) {
            free(nbDistroAA); free(nbDistroAB); free(nbDistroBA); free(nbDistroBB);
        }
    }

    if (doBondedCen==1) {
        free(n_distro_bonded_to_cen_9B); free(n_distro_bonded_to_cen_9K);
        free(n_distro_bonded_to_cen_10B); free(n_distro_bonded_to_cen_10K); free(n_distro_bonded_to_cen_10W);
        free(n_distro_bonded_to_cen_11A); free(n_distro_bonded_to_cen_11B); free(n_distro_bonded_to_cen_11C); free(n_distro_bonded_to_cen_11W);
        free(n_distro_bonded_to_cen_12A); free(n_distro_bonded_to_cen_12B); free(n_distro_bonded_to_cen_12K);
        free(n_distro_bonded_to_cen_13A); free(n_distro_bonded_to_cen_13B); free(n_distro_bonded_to_cen_13K);
        free(n_distro_bonded_to_cen_FCC); free(n_distro_bonded_to_cen_HCP); free(n_distro_bonded_to_cen_BCC_9); free(n_distro_bonded_to_cen_BCC_15);
    }
    
    if (doClusBLDeviation==1) {
        free(bl_mom_sp3); free(bl_mom_sp3a); free(bl_mom_sp3b); free(bl_mom_sp3c);
        free(bl_mom_sp4); free(bl_mom_sp4a); free(bl_mom_sp4b); free(bl_mom_sp4c);
        free(bl_mom_sp5); free(bl_mom_sp5a); free(bl_mom_sp5b); free(bl_mom_sp5c);
        free(bl_mom_6A); free(bl_mom_6Z); free(bl_mom_7K);
        free(bl_mom_8A); free(bl_mom_8B); free(bl_mom_8K);
        free(bl_mom_9A); free(bl_mom_9B); free(bl_mom_9K);
        free(bl_mom_10A); free(bl_mom_10B); free(bl_mom_10K); free(bl_mom_10W);
        free(bl_mom_11A); free(bl_mom_11B); free(bl_mom_11C); free(bl_mom_11E); free(bl_mom_11F); free(bl_mom_11W);
        free(bl_mom_12A); free(bl_mom_12B); free(bl_mom_12D); free(bl_mom_12E); free(bl_mom_12K);
        free(bl_mom_13A); free(bl_mom_13B); free(bl_mom_13K);
        free(bl_mom_FCC); free(bl_mom_HCP); free(bl_mom_BCC_9); free(bl_mom_BCC_15);
    }

    if (doClusBLDistros==1) {
        free(BLDistrosp3); free(BLDistrosp3a); free(BLDistrosp3b); free(BLDistrosp3c);
        free(BLDistrosp4); free(BLDistrosp4a); free(BLDistrosp4b); free(BLDistrosp4c);
        free(BLDistrosp5); free(BLDistrosp5a); free(BLDistrosp5b); free(BLDistrosp5c);
        free(BLDistro6A); free(BLDistro6Z); free(BLDistro7K);
        free(BLDistro8A); free(BLDistro8B); free(BLDistro8K);
        free(BLDistro9A); free(BLDistro9B); free(BLDistro9K);
        free(BLDistro10A); free(BLDistro10B); free(BLDistro10K); free(BLDistro10W);
        free(BLDistro11A); free(BLDistro11B); free(BLDistro11C); free(BLDistro11E); free(BLDistro11F); free(BLDistro11W);
        free(BLDistro12A); free(BLDistro12B); free(BLDistro12D); free(BLDistro12E); free(BLDistro12K);
        free(BLDistro13A); free(BLDistro13B); free(BLDistro13K);
        free(BLDistroFCC); free(BLDistroHCP); free(BLDistroBCC_9); free(BLDistroBCC_15);
    }
}

void Setup_InitPotentialVars(char *filename) { // Initialize variables for potential energy calculations
    char errMsg[1000], input[10000];
    int i, j, throwaway;
    double sigma;
    int ix, iy, iz, imap;
    FILE *fin;
    FILE *fSigma;
    
    printf("d%d reading potential parameters from %s\n",rank,filename);
    fin=fopen(filename,"r");
    if (fin==NULL)  {
        sprintf(errMsg,"Setup_InitPotentialVars() : Error opening file %s",filename);   // Always test file open
        Error(errMsg);
    }
    
    part_pot = malloc(N*sizeof(double));    if (part_pot==NULL) { sprintf(errMsg,"Setup_InitPotentialVars(): part_pot[] malloc out of memory\n");   Error(errMsg); }
    for (i=0; i<N; i++) part_pot[i]=0.0;
    
    sigma_AB=Setup_GetFirstDoubleFromLine(fin);
    sigma_BB=Setup_GetFirstDoubleFromLine(fin);
    epsilon_AB=Setup_GetFirstDoubleFromLine(fin);
    epsilon_BB=Setup_GetFirstDoubleFromLine(fin);
    rcut=Setup_GetFirstDoubleFromLine(fin);
    rcut_AB=Setup_GetFirstDoubleFromLine(fin);
    rcut_BB=Setup_GetFirstDoubleFromLine(fin);
    
    if (fgets(input,10000,fin)==NULL) {
        sprintf(errMsg,"Setup_InitPotentialVars(): empty %s file!",filename);   // Always test file open
        Error(errMsg);
    }
    sprintf(fSigmaName,strtok(input," "));
    
    RHO0=Setup_GetFirstDoubleFromLine(fin);
    mepsilon=Setup_GetFirstDoubleFromLine(fin);
    mcut=Setup_GetFirstDoubleFromLine(fin);
    
    KAPPA=Setup_GetFirstDoubleFromLine(fin);
    yukepsilon=Setup_GetFirstDoubleFromLine(fin);
    yukcut=Setup_GetFirstDoubleFromLine(fin);
    
    ipl_exp=Setup_GetFirstDoubleFromLine(fin);
    ipl_pre=Setup_GetFirstDoubleFromLine(fin);
    
    cubic_a_AA=Setup_GetFirstDoubleFromLine(fin);
    cubic_a_AB=Setup_GetFirstDoubleFromLine(fin);
    cubic_a_BB=Setup_GetFirstDoubleFromLine(fin);
    
    TWOPARTRHO=Setup_GetFirstDoubleFromLine(fin);
    NOPOTENTIALBINS=Setup_GetFirstDoubleFromLine(fin);
    WRITEPOTMAX=Setup_GetFirstDoubleFromLine(fin);
    WRITEPOTMIN=Setup_GetFirstDoubleFromLine(fin);
    
    printf("\n");
    if (WHICHPOTENTIAL==0) {
        printf("d%d %d: (Binary) Lennard-Jones potential\n",rank,WHICHPOTENTIAL);
        printf("d%d LJ simga_AA 1.0 sigma_AB %lg sigma_BB %lg\n",rank,sigma_AB,sigma_BB);
        printf("d%d LJ epsilon_AA 1.0 epsilon_AB %lg epsilon_BB %lg\n",rank,epsilon_AB,epsilon_BB);
        printf("d%d LJ rcut_AA %lg rcut_AB %lg rcut_BB %lg\n",rank,rcut,rcut_AB,rcut_BB);
        
        sigma_AB6=pow(sigma_AB,6.0);
        sigma_BB6=pow(sigma_BB,6.0);
        
        rcut2=rcut*rcut;
        rcut_AB2=rcut_AB*rcut_AB;
        rcut_BB2=rcut_BB*rcut_BB;
        
        quarteruTail_AA = pow(rcut, -12.0) - pow(rcut, -6.0);   // tail on lennard jones potential because of cut off
        quarteruTail_AB = epsilon_AB*(sigma_AB6*sigma_AB6*pow(rcut_AB, -12.0) - sigma_AB6*pow(rcut_AB, -6.0));  // tail on lennard jones potential because of cut off
        quarteruTail_BB = epsilon_BB*(sigma_BB6*sigma_BB6*pow(rcut_BB, -12.0) - sigma_BB6*pow(rcut_BB, -6.0));  // tail on lennard jones potential because of cut off
    }
    else if (WHICHPOTENTIAL==1) {
        printf("d%d %d: Stoddard-Ford cut-off (Binary) Lennard-Jones potential\n",rank,WHICHPOTENTIAL);
        printf("d%d SFLJ simga_AA 1.0 sigma_AB %lg sigma_BB %lg\n",rank,sigma_AB,sigma_BB);
        printf("d%d SFLJ epsilon_AA 1.0 epsilon_AB %lg epsilon_BB %lg\n",rank,epsilon_AB,epsilon_BB);
        printf("d%d SFLJ rcut_AA %lg rcut_AB %lg rcut_BB %lg\n",rank,rcut,rcut_AB,rcut_BB);
        
        sigma_AB6=pow(sigma_AB,6.0);
        sigma_BB6=pow(sigma_BB,6.0);
        
        rcut2=rcut*rcut;
        rcut_AB2=rcut_AB*rcut_AB;
        rcut_BB2=rcut_BB*rcut_BB;
        
        quarteruTail_AA = pow(rcut, -12.0) - pow(rcut, -6.0);   // tail on lennard jones potential because of cut off
        quarteruTail_AB = epsilon_AB*(sigma_AB6*sigma_AB6*pow(rcut_AB, -12.0) - sigma_AB6*pow(rcut_AB, -6.0));  // tail on lennard jones potential because of cut off
        quarteruTail_BB = epsilon_BB*(sigma_BB6*sigma_BB6*pow(rcut_BB, -12.0) - sigma_BB6*pow(rcut_BB, -6.0));  // tail on lennard jones potential because of cut off
        stoddardford_AA = (6.0*pow(rcut, -12.0) - 3.0*pow(rcut, -6.0))/(rcut2);
        stoddardford_AB = epsilon_AB*(6.0*sigma_AB6*sigma_AB6*pow(rcut_AB, -12.0) - 3.0*sigma_AB6*pow(rcut_AB, -6.0))/(rcut_AB2);
        stoddardford_BB = epsilon_BB*(6.0*sigma_BB6*sigma_BB6*pow(rcut_BB, -12.0) - 3.0*sigma_BB6*pow(rcut_BB, -6.0))/(rcut_BB2);
    }
    else if (WHICHPOTENTIAL==2) {
        printf("d%d %d: Morse+Yukawa potential\n",rank,WHICHPOTENTIAL);
        printf("d%d Morse rho0 %lg\n",rank,RHO0);
        printf("d%d Morse mepsilon %lg\n",rank,mepsilon);
        printf("d%d Morse mcut %lg\n",rank,mcut);
        
        printf("d%d Yukawa kappa %lg\n",rank,KAPPA);
        printf("d%d Yukawa yukepsilon %lg\n",rank,yukepsilon);
        printf("d%d Yukawa yukcut %lg\n",rank,yukcut);
        
        psigma = malloc(N*sizeof(double));  if (psigma==NULL) { sprintf(errMsg,"Setup_InitPotentialVars(): psigma[] malloc out of memory\n");   Error(errMsg); }
        for (i=0; i<N; i++) psigma[i]=0.0;
        
        printf("d%d reading particle diameters from %s\n",rank,fSigmaName);
        fSigma=fopen(fSigmaName,"r");
        if (fSigma==NULL)  {
            sprintf(errMsg,"Setup_InitPotentialVars() : Error opening file %s",fSigmaName); // Always test file open
            Error(errMsg);
        }
        
        for(i=0; i<N; ++i) {
            fscanf(fSigma,"%d   %lg",&throwaway,&sigma);
            psigma[i] = sigma;
        }
        
        fclose(fSigma);
        
        mcut2=mcut*mcut;
        yukcut2=yukcut*yukcut;
        
        uMorseTail = (mepsilon * exp(RHO0*(1.0-mcut)) * (exp(RHO0*(1.0-mcut)) - 2.0) ); // tail on morse potential because of cut off
        uYukTail = (yukepsilon * exp(-KAPPA*(yukcut-1.0)) / yukcut);    // tail on yukawa potential because of cut off
        if (yukcut==0.0) uYukTail=0.0;
    }
    else if (WHICHPOTENTIAL==3) {
        printf("d%d %d: Not currently implemented potential\n",rank,WHICHPOTENTIAL);
    }
    else if (WHICHPOTENTIAL==4) {
        printf("d%d %d: Inverse-power law potential\n",rank,WHICHPOTENTIAL);
        printf("d%d IPL exponent %lg\n",rank,ipl_exp);
        printf("d%d IPL epsilon_AA 1.0 epsilon_AB %lg epsilon_BB %lg\n",rank,epsilon_AB,epsilon_BB);
        printf("d%d IPL prefactor %lg\n",rank,ipl_pre);
        
        half_ipl_exp=ipl_exp/2.0;
        
        sigma_AB2=pow(sigma_AB,2.0);
        sigma_BB2=pow(sigma_BB,2.0);
        
        rcut2=rcut*rcut;
        rcut_AB2=rcut_AB*rcut_AB;
        rcut_BB2=rcut_BB*rcut_BB;
        
        uTail_AA = pow((1.0/rcut), ipl_exp);    // tail on lennard jones potential because of cut off
        uTail_AB = epsilon_AB*pow((sigma_AB/rcut_AB), ipl_exp); // tail on lennard jones potential because of cut off
        uTail_BB = epsilon_BB*pow((sigma_BB/rcut_BB), ipl_exp); // tail on lennard jones potential because of cut off
    }
    else if (WHICHPOTENTIAL==5) {
        sigma_AB6=pow(sigma_AB,6.0);
        sigma_BB6=pow(sigma_BB,6.0);
        
        rcut2=rcut*rcut;
        rcut_AB2=rcut_AB*rcut_AB;
        rcut_BB2=rcut_BB*rcut_BB;
        
        cubic_a_AA2=cubic_a_AA*cubic_a_AA;
        cubic_a_AB2=cubic_a_AB*cubic_a_AB;
        cubic_a_BB2=cubic_a_BB*cubic_a_BB;
        
        cubicB_AA=-24.0*(2.0*pow(cubic_a_AA, -13.0)-pow(cubic_a_AA, -7.0))/(-3.0*pow((rcut-cubic_a_AA),2.0));
        cubicB_AB=-24.0*epsilon_AB*(2.0*sigma_AB6*sigma_AB6*pow(cubic_a_AB, -13.0)-sigma_AB6*pow(cubic_a_AB, -7.0))/(-3.0*pow((rcut_AB-cubic_a_AB),2.0));
        cubicB_BB=-24.0*epsilon_BB*(2.0*sigma_BB6*sigma_BB6*pow(cubic_a_BB, -13.0)-sigma_BB6*pow(cubic_a_BB, -7.0))/(-3.0*pow((rcut_BB-cubic_a_BB),2.0));
        
        cubicA_AA=-4.0*(pow(cubic_a_AA, -12.0)-pow(cubic_a_AA, -6.0))+cubicB_AA*pow((rcut-cubic_a_AA),3.0);
        cubicA_AB=-4.0*epsilon_AB*(sigma_AB6*sigma_AB6*pow(cubic_a_AB, -12.0)-sigma_AB6*pow(cubic_a_AB, -6.0))+cubicB_AB*pow((rcut_AB-cubic_a_AB),3.0);
        cubicA_BB=-4.0*epsilon_BB*(sigma_BB6*sigma_BB6*pow(cubic_a_BB, -12.0)-sigma_BB6*pow(cubic_a_BB, -6.0))+cubicB_BB*pow((rcut_BB-cubic_a_BB),3.0);
        
        printf("d%d Binary Potential\n",rank);
        printf("d%d sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",rank,sigma_AB,sigma_BB);
        printf("d%d epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",rank,epsilon_AB,epsilon_BB);
        printf("d%d N %d N_A %d N_B %d\n",rank,N,NA,N-NA);
        
        printf("d%d cubic spline Truncation and Smoothing\n",rank);
        printf("d%d Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rank,rcut,rcut_AB,rcut_BB);
        printf("d%d LJ to cubic cross-over lengths: cubic_a_AA %.15lg cubic_a_AB %.15lg cubic_a_BB %.15lg\n",rank,cubic_a_AA,cubic_a_AB,cubic_a_BB);
        printf("d%d LJ to cubic A params: cubicA_AA %.15lg cubicA_AB %.15lg cubicA_BB %.15lg\n",rank,cubicA_AA,cubicA_AB,cubicA_BB);
        printf("d%d LJ to cubic B params: cubicB_AA %.15lg cubicB_AB %.15lg cubicB_BB %.15lg\n",rank,cubicB_AA,cubicB_AB,cubicB_BB);
    }
    else if (WHICHPOTENTIAL==6) {
        printf("d%d %d: Stoddard-Ford cut-off (Binary) Inverse-power law potential\n",rank,WHICHPOTENTIAL);
        printf("d%d SFIPL exponent %lg\n",rank,ipl_exp);
        printf("d%d SFIPL epsilon_AA 1.0 epsilon_AB %lg epsilon_BB %lg\n",rank,epsilon_AB,epsilon_BB);
        printf("d%d SFIPL prefactor %lg\n",rank,ipl_pre);
        printf("d%d SFIPL rcut_AA %lg rcut_AB %lg rcut_BB %lg\n",rank,rcut,rcut_AB,rcut_BB);
        
        half_ipl_exp=ipl_exp/2.0;
        
        sigma_AB2=pow(sigma_AB,2.0);
        sigma_BB2=pow(sigma_BB,2.0);
        
        rcut2=rcut*rcut;
        rcut_AB2=rcut_AB*rcut_AB;
        rcut_BB2=rcut_BB*rcut_BB;
        
        uTail_AA = pow((1.0/rcut), ipl_exp);    // tail on lennard jones potential because of cut off
        uTail_AB = epsilon_AB*pow((sigma_AB/rcut_AB), ipl_exp); // tail on lennard jones potential because of cut off
        uTail_BB = epsilon_BB*pow((sigma_BB/rcut_BB), ipl_exp); // tail on lennard jones potential because of cut off
        
        stoddardford_AA = (half_ipl_exp*pow(rcut, -ipl_exp))/(rcut2);
        stoddardford_AB = epsilon_AB*(half_ipl_exp*pow((sigma_AB/rcut_AB), ipl_exp))/(rcut_AB2);
        stoddardford_BB = epsilon_BB*(half_ipl_exp*pow((sigma_BB/rcut_BB), ipl_exp))/(rcut_BB2);    
    }
    else if (WHICHPOTENTIAL==7) {
        printf("d%d %d: CRVT potential with min at r=1 and well depth \varepsilon=1\n",rank,WHICHPOTENTIAL);
        rcut2=rcut*rcut;
    }
    fclose(fin);
    
    if (USELIST==1) {
        M_pot = (int)(side/rcut);   // number of cells along box side
        if (M_pot<3) Error("Setup_InitPotentialVars(): M_pot<3, too few cells");
        ncells_pot = M_pot*M_pot*M_pot; // total number of cells
        
        map_pot=malloc(((13*ncells_pot)+1)*sizeof(int));    if (map_pot==NULL) { sprintf(errMsg,"Setup_InitPotentialVars(): map_pot[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<(13*ncells_pot+1); j++) map_pot[j]=0;
        head_pot=malloc((ncells_pot+1)*sizeof(int));    if (head_pot==NULL) { sprintf(errMsg,"Setup_InitPotentialVars(): head_pot[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<(ncells_pot+1); j++) head_pot[j]=0;
        llist_pot=malloc((N+1)*sizeof(int));    if (llist_pot==NULL) { sprintf(errMsg,"Setup_InitPotentialVars(): llist_pot[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<(N+1); j++) llist_pot[j]=0;
        
        cellSide_pot = side/M_pot;  // length of cells
        invcellSide_pot = 1.0/cellSide_pot; // invcellSide
        
        printf("d%d potential cell list m_pot %d ncells_pot %d cellside_pot %.15lg\n", rank,M_pot, ncells_pot, cellSide_pot);
        // routine to create the thirteen nearest neighbours array map[] of each cell 
        for (iz=1; iz<=M_pot; iz++) {
            for (iy=1; iy<=M_pot; iy++) {
                for (ix=1; ix<=M_pot; ix++) {
                    imap = (icell_pot(ix,iy,iz)-1)*13;
                    map_pot[imap+1 ]=icell_pot(ix+1,iy  ,iz );
                    map_pot[imap+2 ]=icell_pot(ix+1,iy+1,iz );
                    map_pot[imap+3 ]=icell_pot(ix    ,iy+1,iz   );
                    map_pot[imap+4 ]=icell_pot(ix-1 ,iy+1,iz    );
                    map_pot[imap+5 ]=icell_pot(ix+1,iy  ,iz-1   );
                    map_pot[imap+6 ]=icell_pot(ix+1,iy+1,iz-1   );
                    map_pot[imap+7 ]=icell_pot(ix    ,iy+1,iz-1 );
                    map_pot[imap+8 ]=icell_pot(ix-1 ,iy+1,iz-1  );
                    map_pot[imap+9 ]=icell_pot(ix+1,iy  ,iz+1   );
                    map_pot[imap+10]=icell_pot(ix+1,iy+1,iz+1   );
                    map_pot[imap+11]=icell_pot(ix    ,iy+1,iz+1 );
                    map_pot[imap+12]=icell_pot(ix-1 ,iy+1,iz+1);
                    map_pot[imap+13]=icell_pot(ix    ,iy    ,iz+1   );
                }
            }
        }
        
        free(map_pot);
        free(head_pot);
        free(llist_pot);
    }
    
    av_potential=av_pot_check=0.0;
    
    av_pot_sp3=av_pot_sp3a=av_pot_sp3b=av_pot_sp3c=0.0;
    av_pot_sp4=av_pot_sp4a=av_pot_sp4b=av_pot_6A=0.0;
    av_pot_sp5=av_pot_sp5a=av_pot_sp5b=av_pot_sp5c=0.0;
    av_pot_6Z=av_pot_7K=0.0;
    av_pot_8A=av_pot_8B=av_pot_8K=0.0;
    av_pot_9A=av_pot_9B=av_pot_9K=0.0;
    av_pot_10A=av_pot_10B=av_pot_10K=av_pot_10W=0.0;
    av_pot_11A=av_pot_11B=av_pot_11C=av_pot_11E=av_pot_11F=av_pot_11W=0.0;
    av_pot_12A=av_pot_12B=av_pot_12D=av_pot_12E=av_pot_12K=0.0;
    av_pot_13A=av_pot_13B=av_pot_13K=0.0;
    av_pot_FCC=av_pot_HCP=av_pot_BCC_9=av_pot_BCC_15=0.0;

    av_pot_cen_9B=av_pot_cen_9K=0.0;
    av_pot_cen_10B=av_pot_cen_10K=av_pot_cen_10W=0.0;
    av_pot_cen_11A=av_pot_cen_11B=av_pot_cen_11C=av_pot_cen_11W=0.0;
    av_pot_cen_12A=av_pot_cen_12B=av_pot_cen_12K=0.0;
    av_pot_cen_13A=av_pot_cen_13B=av_pot_cen_13K=0.0;
    av_pot_cen_FCC=av_pot_cen_HCP=av_pot_cen_BCC_9=av_pot_cen_BCC_15=0.0;
    
    cnt_av_pot_cen_9B=cnt_av_pot_cen_9K=0;
    cnt_av_pot_cen_10B=cnt_av_pot_cen_10K=cnt_av_pot_cen_10W=0;
    cnt_av_pot_cen_11A=cnt_av_pot_cen_11B=cnt_av_pot_cen_11C=cnt_av_pot_cen_11W=0;
    cnt_av_pot_cen_12A=cnt_av_pot_cen_12B=cnt_av_pot_cen_12K=0;
    cnt_av_pot_cen_13A=cnt_av_pot_cen_13B=cnt_av_pot_cen_13K=0;
    cnt_av_pot_cen_FCC=cnt_av_pot_cen_HCP=cnt_av_pot_cen_BCC_9=cnt_av_pot_cen_BCC_15=0;
    
    av_pot_shell_9B=av_pot_shell_9K=0.0;
    av_pot_shell_10B=av_pot_shell_10K=av_pot_shell_10W=0.0;
    av_pot_shell_11A=av_pot_shell_11B=av_pot_shell_11C=av_pot_shell_11W=0.0;
    av_pot_shell_12A=av_pot_shell_12B=av_pot_shell_12K=0.0;
    av_pot_shell_13A=av_pot_shell_13B=av_pot_shell_13K=0.0;
    av_pot_shell_FCC=av_pot_shell_HCP=av_pot_shell_BCC_9=av_pot_shell_BCC_15=0.0;
    
    cnt_av_pot_shell_9B=cnt_av_pot_shell_9K=0;
    cnt_av_pot_shell_10B=cnt_av_pot_shell_10K=cnt_av_pot_shell_10W=0;
    cnt_av_pot_shell_11A=cnt_av_pot_shell_11B=cnt_av_pot_shell_11C=cnt_av_pot_shell_11W=0;
    cnt_av_pot_shell_12A=cnt_av_pot_shell_12B=cnt_av_pot_shell_12K=0;
    cnt_av_pot_shell_13A=cnt_av_pot_shell_13B=cnt_av_pot_shell_13K=0;
    cnt_av_pot_shell_FCC=cnt_av_pot_shell_HCP=cnt_av_pot_shell_BCC_9=cnt_av_pot_shell_BCC_15=0;
}

void Setup_FreePotentialVars() {
    
    free(part_pot);
    if (WHICHPOTENTIAL==2) free(psigma);
    if (USELIST==1) {
        free(map_pot);
        free(llist_pot);
        free(head_pot);
    }
}
void Setup_print_U_r() {    // prints details of potential used
    int i, j, storeN, storeRtype0, storeRtype1, storeM_pot, storencells_pot;
    char errMsg[1000], output[1000];
    double rsep, storeSide, storeHalfSide, storecellSide_pot, storeinvcellSide_pot;
    int ix, iy, iz;
    int imap;
    FILE *potForceFile;
    
    storeM_pot=storencells_pot=0;
    storecellSide_pot=storeinvcellSide_pot=0.0; 
    storeN=N;
    N=2;
    storeSide=side;
    side = pow(N/TWOPARTRHO,1.0/3.0);
    storeHalfSide=halfSide;
    halfSide=side/2.0;
    if (USELIST==1) {
        storeM_pot=M_pot;
        M_pot = (int)(side/rcut);   // number of cells along box side
        if (M_pot<3) Error("Setup_print_U_r(): M_pot<3, too few cells");
        storencells_pot = ncells_pot;
        ncells_pot = M_pot*M_pot*M_pot; // total number of cells
        storecellSide_pot=cellSide_pot;
        cellSide_pot = side/M_pot;  // length of cells
        storeinvcellSide_pot=invcellSide_pot;
        invcellSide_pot = 1.0/cellSide_pot; // invcellSide
        map_pot=malloc(((13*ncells_pot)+1)*sizeof(int));    if (map_pot==NULL) { sprintf(errMsg,"Setup_print_U_r(): map_pot[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<(13*ncells_pot+1); j++) map_pot[j]=0;
        head_pot=malloc((ncells_pot+1)*sizeof(int));    if (head_pot==NULL) { sprintf(errMsg,"Setup_print_U_r(): head_pot[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<(ncells_pot+1); j++) head_pot[j]=0;
        llist_pot=malloc((N+1)*sizeof(int));    if (llist_pot==NULL) { sprintf(errMsg,"Setup_print_U_r(): llist_pot[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<(N+1); j++) llist_pot[j]=0;
        
        //printf("d%d m %d ncells %d cellside %.15lg\n", rank,M, ncells, cellSide);
        // routine to create the thirteen nearest neighbours array map[] of each cell 
        for (iz=1; iz<=M_pot; iz++) {
            for (iy=1; iy<=M_pot; iy++) {
                for (ix=1; ix<=M_pot; ix++) {
                    imap = (icell_pot(ix,iy,iz)-1)*13;
                    map_pot[imap+1 ]=icell_pot(ix+1,iy  ,iz );
                    map_pot[imap+2 ]=icell_pot(ix+1,iy+1,iz );
                    map_pot[imap+3 ]=icell_pot(ix    ,iy+1,iz   );
                    map_pot[imap+4 ]=icell_pot(ix-1 ,iy+1,iz    );
                    map_pot[imap+5 ]=icell_pot(ix+1,iy  ,iz-1   );
                    map_pot[imap+6 ]=icell_pot(ix+1,iy+1,iz-1   );
                    map_pot[imap+7 ]=icell_pot(ix    ,iy+1,iz-1 );
                    map_pot[imap+8 ]=icell_pot(ix-1 ,iy+1,iz-1  );
                    map_pot[imap+9 ]=icell_pot(ix+1,iy  ,iz+1   );
                    map_pot[imap+10]=icell_pot(ix+1,iy+1,iz+1   );
                    map_pot[imap+11]=icell_pot(ix    ,iy+1,iz+1 );
                    map_pot[imap+12]=icell_pot(ix-1 ,iy+1,iz+1);
                    map_pot[imap+13]=icell_pot(ix    ,iy    ,iz+1   );
                }
            }
        }
    }
    
    storeRtype0=rtype[0];
    storeRtype1=rtype[1];
    
    rsep=(WRITEPOTMAX-WRITEPOTMIN)/NOPOTENTIALBINS;
    x[0]=y[0]=z[0]=x[1]=y[1]=0.0;
    
    if (WHICHPOTENTIAL==0) {
        if (USELIST==0) sprintf(output,"d%d_tcc.BLJ.Ur",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.BLJ.Ur",rank);
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
        fprintf(potForceFile,"%s Binary Lennard-Jones type potenial\n",output);
        fprintf(potForceFile,"r U_AA(r) U_AB(r) U_BB(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
            
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) BLJ();
            else {
                links_pot();
                listBLJ();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=1;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJ();
            else {
                links_pot();
                listBLJ();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJ();
            else {
                links_pot();
                listBLJ();
            }
            fprintf(potForceFile,"  %lg\n",potential);
        }
        
        fprintf(potForceFile,"\nBinary Potential\n");
        fprintf(potForceFile,"sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",sigma_AB,sigma_BB);
        fprintf(potForceFile,"epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",epsilon_AB,epsilon_BB);
        fprintf(potForceFile,"Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rcut,rcut_AB,rcut_BB);
        fprintf(potForceFile,"Pot shifts: uTail_AA %.15lg uTail_AB %.15lg uTail_BB %.15lg\n",4.0*quarteruTail_AA,4.0*quarteruTail_AB,4.0*quarteruTail_BB);
    }
    
    else if (WHICHPOTENTIAL==1) {
        if (USELIST==0) sprintf(output,"d%d_tcc.SF.BLJ.Ur",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.SF.BLJ.Ur",rank);
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
        fprintf(potForceFile,"%s Stoddard-ford Binary Lennard-Jones type potenial\n",output);
        fprintf(potForceFile,"r U_AA(r) U_AB(r) U_BB(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
            
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) BLJSF();
            else {
                links_pot();
                listBLJSF();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=1;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJSF();
            else {
                links_pot();
                listBLJSF();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJSF();
            else {
                links_pot();
                listBLJSF();
            }
            fprintf(potForceFile,"  %lg\n",potential);
        }
        
        fprintf(potForceFile,"\nBinary Potential\n");
        fprintf(potForceFile,"sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",sigma_AB,sigma_BB);
        fprintf(potForceFile,"epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",epsilon_AB,epsilon_BB);
        fprintf(potForceFile,"Stoddard-Ford Truncation and Smoothing\n");
        fprintf(potForceFile,"Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rcut,rcut_AB,rcut_BB);
        fprintf(potForceFile,"Pot shifts: SFuTail_AA %.15lg SFuTail_AB %.15lg SFuTail_BB %.15lg\n",4.0*quarteruTail_AA,4.0*quarteruTail_AB,4.0*quarteruTail_BB);
        fprintf(potForceFile,"Force shifts: SFfTail_AA %.15lg SFfTail_AB %.15lg SFfTail_BB %.15lg\n",8.0*stoddardford_AA,8.0*stoddardford_AB,8.0*stoddardford_BB);
    }
    
    else if (WHICHPOTENTIAL==2) {
        if (USELIST==0) sprintf(output,"d%d_tcc.Mor.Yuk.Ur",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.Mor.Yuk.Ur",rank);
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
        fprintf(potForceFile,"%s Morse+Yukawa potenials\n",output);
        fprintf(potForceFile,"r U(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
            
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if(USELIST==0) MorYuk();
            else {
                links_pot();
                listMorYuk();
            }
            fprintf(potForceFile,"  %lg\n",potential);
        }
        fprintf(potForceFile,"\nMorse: BEM %.15lg RHO0 %.15lg\n",mepsilon,RHO0);
        fprintf(potForceFile,"Yuk: BEY %.15lg Kappa %.15lg \n",yukepsilon,KAPPA);
        fprintf(potForceFile,"Cut-offs: Morse %.15lg Yukawa %.15lg \n",mcut,yukcut);
        fprintf(potForceFile,"Pot shifts: uMorseTail %.15lg uYukTail %.15lg\n",uMorseTail,uYukTail);
    }
    
    else if (WHICHPOTENTIAL==4) {
        if (USELIST==0) sprintf(output,"d%d_tcc.BIPL.Ur",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.BIPL.Ur",rank);
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
        fprintf(potForceFile,"%s Binary IPL type potenial\n",output);
        fprintf(potForceFile,"r U_AA(r) U_AB(r) U_BB(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
            
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) BIPL();
            else {
                links_pot();
                listBIPL();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=1;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BIPL();
            else {
                links_pot();
                listBIPL();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BIPL();
            else {
                links_pot();
                listBIPL();
            }
            fprintf(potForceFile,"  %lg\n",potential);
        }
        
        fprintf(potForceFile,"\nBinary IPL Potential\n");
        fprintf(potForceFile,"ipl_prefactor %.15lg ipl_exponent %.15lg half_ipl_exponent %.15lg\n",ipl_pre,ipl_exp,half_ipl_exp);
        fprintf(potForceFile,"sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",sigma_AB,sigma_BB);
        fprintf(potForceFile,"epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",epsilon_AB,epsilon_BB);
        fprintf(potForceFile,"Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rcut,rcut_AB,rcut_BB);
        fprintf(potForceFile,"Pot shifts: uTail_AA %.15lg uTail_AB %.15lg uTail_BB %.15lg\n",ipl_pre*uTail_AA,ipl_pre*uTail_AB,ipl_pre*uTail_BB);
    }
    
    else if (WHICHPOTENTIAL==5) {
        if (USELIST==0) sprintf(output,"d%d_tcc.BLJ_WCA_s.UrFr",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.BLJ_WCA_s.UrFr",rank);
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
        fprintf(potForceFile,"%s Binary Lennard-Jones WCA cubic splined type potenial\n",output);
        fprintf(potForceFile,"r U_AA(r) U_AB(r) U_BB(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
            
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) BLJ_WCA_s();
            else {
                links_pot();
                listBLJ_WCA_s();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=1;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJ_WCA_s();
            else {
                links_pot();
                listBLJ_WCA_s();
            }
            fprintf(potForceFile,"  %lg",potential);
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) BLJ_WCA_s();
            else {
                links_pot();
                listBLJ_WCA_s();
            }
            fprintf(potForceFile,"  %lg\n",potential);
        }
        
        fprintf(potForceFile,"\nd%d Binary Lennard-Jones WCA cubic splined Potential\n",rank);
        fprintf(potForceFile,"d%d sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",rank,sigma_AB,sigma_BB);
        fprintf(potForceFile,"d%d epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",rank,epsilon_AB,epsilon_BB);
        fprintf(potForceFile,"d%d cubic spline Truncation and Smoothing\n",rank);
        fprintf(potForceFile,"d%d Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rank,rcut,rcut_AB,rcut_BB);
        fprintf(potForceFile,"d%d LJ to cubic cross-over lengths: cubic_a_AA %.15lg cubic_a_AB %.15lg cubic_a_BB %.15lg\n",rank,cubic_a_AA,cubic_a_AB,cubic_a_BB);
        fprintf(potForceFile,"d%d LJ to cubic A params: cubicA_AA %.15lg cubicA_AB %.15lg cubicA_BB %.15lg\n",rank,cubicA_AA,cubicA_AB,cubicA_BB);
        fprintf(potForceFile,"d%d LJ to cubic B params: cubicB_AA %.15lg cubicB_AB %.15lg cubicB_BB %.15lg\n",rank,cubicB_AA,cubicB_AB,cubicB_BB);
    }
    
    else if (WHICHPOTENTIAL==6) {
        if (USELIST==0) sprintf(output,"d%d_tcc.SFBIPL.UrFr",rank);
        else sprintf(output,"d%d_tcc.CELL_LIST.SFBIPL.UrFr",rank);
        potForceFile=fopen(output,"w");
        fprintf(potForceFile,"%s Stoddard-Ford Binary IPL type potenial\n",output);
        fprintf(potForceFile,"r U_AA(r) U_AB(r) U_BB(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
        
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) SFBIPL();
            else if (USELIST==1) {
                links_pot();
                listSFBIPL();
            }
            fprintf(potForceFile,"  %.15lg",potential);
            rtype[0]=1;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) SFBIPL();
            else if (USELIST==1) {
                links_pot();
                listSFBIPL();
            }
            fprintf(potForceFile,"  %.15lg",potential);
            rtype[0]=2;
            rtype[1]=2;
            potential=0.0;
            if (USELIST==0) SFBIPL();
            else if (USELIST==1) {
                links_pot();
                listSFBIPL();
            }
            fprintf(potForceFile,"  %.15lg\n",potential);
        }
            
        fprintf(potForceFile,"\nd%d Stoddard-Ford Binary IPL Potential\n",rank);
        fprintf(potForceFile,"d%d ipl_prefactor %.15lg ipl_exponent %.15lg half_ipl_exponent %.15lg\n",rank,ipl_pre,ipl_exp,half_ipl_exp);
        fprintf(potForceFile,"d%d sigma_AA 1.0 sigma_AB %.15lg sigma_BB %.15lg\n",rank,sigma_AB,sigma_BB);
        fprintf(potForceFile,"d%d epsilon_AA 1.0 epsilon_AB %.15lg epsilon_BB %.15lg\n",rank,epsilon_AB,epsilon_BB);
        fprintf(potForceFile,"d%d Stoddard-Ford Truncation and Smoothing\n",rank);
        fprintf(potForceFile,"d%d Cut-offs: rcut_AA %.15lg rcut_AB %.15lg rcut_BB %.15lg\n",rank,rcut,rcut_AB,rcut_BB);
        fprintf(potForceFile,"d%d Pot shifts: SFuTail_AA %.15lg SFuTail_AB %.15lg SFuTail_BB %.15lg\n",rank,ipl_pre*uTail_AA,ipl_pre*uTail_AB,ipl_pre*uTail_BB);
        fprintf(potForceFile,"d%d Force shifts: SFfTail_AA %.15lg SFfTail_AB %.15lg SFfTail_BB %.15lg\n",rank,2.0*ipl_pre*rcut*stoddardford_AA,2.0*ipl_pre*rcut_AB*stoddardford_AB,2.0*ipl_pre*rcut_BB*stoddardford_BB);
    }
    
    else if (WHICHPOTENTIAL==7) {
        if (USELIST==0) sprintf(output,"d%d_tcc.CRVT.UrFr",rank);
        //else sprintf(output,"d%d_tcc.CELL_LIST.SFBIPL.UrFr",rank);
        potForceFile=fopen(output,"w");
        fprintf(potForceFile,"%s CRVTpotenial\n",output);
        fprintf(potForceFile,"r U(r)\n");
        for (i=0;i<=NOPOTENTIALBINS;i++) {
            z[1]=WRITEPOTMIN+i*rsep;
            fprintf(potForceFile,"%lg",z[1]);
        
            rtype[0]=1;
            rtype[1]=1;
            potential=0.0;
            if (USELIST==0) CRVT();
            /*else if (USELIST==1) {
                links_pot();
                listCRVT();
            }*/
            fprintf(potForceFile,"  %.15lg",potential);
            
        }
            
        fprintf(potForceFile,"\nd%d CRVT Potential\n",rank);
        fprintf(potForceFile,"d%d min r=1 well depth \varepsilon=1\n",rank);
    }
    
    else {
        sprintf(errMsg,"Setup_print_U_r(): WHICHPOTENTIAL %d not recognized\n",WHICHPOTENTIAL);
        sprintf(output,"Setup_print_U_r.dummy");
        potForceFile=fopen(output,"w");
        if (potForceFile==NULL)  {
            sprintf(errMsg,"Setup_print_U_r(): Error opening file %s",output);  // Always test file open
            Error(errMsg);
        }
    }
    fclose(potForceFile);
    printf("d%d Written %s\n",rank,output);
    
    potential=0.0;
    N=storeN;
    side=storeSide;
    halfSide=storeHalfSide;
    
    rtype[0]=storeRtype0;
    rtype[1]=storeRtype1;
    
    if (USELIST==1) {
        M_pot=storeM_pot;
        ncells_pot=storencells_pot;
        cellSide_pot=storecellSide_pot;
        invcellSide_pot=storeinvcellSide_pot;
    }   
}

int icell_pot(int tix, int tiy, int tiz) {  // returns cell number (from 1 to ncells) for given (tix,tiy,tiz) coordinate
    return 1 + (tix-1+M_pot)%M_pot + M_pot*((tiy-1+M_pot)%M_pot) + M_pot*M_pot*((tiz-1+M_pot)%M_pot); 
}

void links_pot() {  // sorts all the particles into cells, result given by head-of-chain and linked list arrays
    int i, ic;
    for (ic=1;ic<=ncells_pot;ic++) head_pot[ic]=0;
    for (i=1;i<=N;i++) {
        ic = 1 + (int)((x[i-1]+ halfSide)*invcellSide_pot) + M_pot*((int)((y[i-1]+halfSide)*invcellSide_pot)) + M_pot*M_pot*((int)((z[i-1]+halfSide)*invcellSide_pot));
        if (ic > ncells_pot || ic <= 0) {
            printf("d%d i %d r_x %lg r_y %lg r_z %lg side %lg halfSide %lg ic %d ncells_pot %d\n",rank,i-1,x[i-1],y[i-1],z[i-1],side,halfSide,ic,ncells_pot);
            Error("links_pot(): ic > ncells_pot i.e. particle coord no longer in simulation box!!\n");
        }
        llist_pot[i]=head_pot[ic];
        head_pot[ic]=i;
    }
}

void Setup_ClusComp() { // zero arrays for cluster compostion analysis
    int i;
    
    for (i=0; i<4; i++) n_distro_sp3[i]=0;
    for (i=0; i<4; i++) n_distro_sp3a[i]=0;
    for (i=0; i<5; i++) n_distro_sp3b[i]=0;
    for (i=0; i<6; i++) n_distro_sp3c[i]=0;
    for (i=0; i<5; i++) n_distro_sp4[i]=0;
    for (i=0; i<5; i++) n_distro_sp4a[i]=0;
    for (i=0; i<6; i++) n_distro_sp4b[i]=0;
    for (i=0; i<7; i++) n_distro_sp4c[i]=0;
    for (i=0; i<7; i++) n_distro_6A[i]=0;
    for (i=0; i<6; i++) n_distro_sp5[i]=0;
    for (i=0; i<6; i++) n_distro_sp5a[i]=0;
    for (i=0; i<7; i++) n_distro_sp5b[i]=0;
    for (i=0; i<8; i++) n_distro_sp5c[i]=0;
    for (i=0; i<7; i++) n_distro_6Z[i]=0;
    for (i=0; i<8; i++) n_distro_7K[i]=0;
    for (i=0; i<9; i++) n_distro_8A[i]=0;
    for (i=0; i<9; i++) n_distro_8B[i]=0;
    for (i=0; i<9; i++) n_distro_8K[i]=0;
    for (i=0; i<10; i++) n_distro_9A[i]=0;
    for (i=0; i<10; i++) n_distro_9B[i]=0;
    for (i=0; i<10; i++) n_distro_9K[i]=0;
    for (i=0; i<11; i++) n_distro_10A[i]=0;
    for (i=0; i<11; i++) n_distro_10B[i]=0;
    for (i=0; i<11; i++) n_distro_10K[i]=0;
    for (i=0; i<11; i++) n_distro_10W[i]=0;
    for (i=0; i<12; i++) n_distro_11A[i]=0;
    for (i=0; i<12; i++) n_distro_11B[i]=0;
    for (i=0; i<12; i++) n_distro_11C[i]=0;
    for (i=0; i<12; i++) n_distro_11E[i]=0;
    for (i=0; i<12; i++) n_distro_11F[i]=0;
    for (i=0; i<12; i++) n_distro_11W[i]=0;
    for (i=0; i<13; i++) n_distro_12A[i]=0;
    for (i=0; i<13; i++) n_distro_12B[i]=0;
    for (i=0; i<13; i++) n_distro_12D[i]=0;
    for (i=0; i<13; i++) n_distro_12E[i]=0;
    for (i=0; i<13; i++) n_distro_12K[i]=0;
    for (i=0; i<14; i++) n_distro_13A[i]=0;
    for (i=0; i<14; i++) n_distro_13B[i]=0;
    for (i=0; i<14; i++) n_distro_13K[i]=0;
    for (i=0; i<14; i++) n_distro_FCC[i]=0;
    for (i=0; i<14; i++) n_distro_HCP[i]=0;
    for (i=0; i<10; i++) n_distro_BCC_9[i]=0;
    for (i=0; i<16; i++) n_distro_BCC_15[i]=0;
    
    for (i=0; i<2; i++) {
        n_distro_cen_9B[i]=0;
        n_distro_cen_9K[i]=0;
        n_distro_cen_10B[i]=0;
        n_distro_cen_10K[i]=0;
        n_distro_cen_10W[i]=0;
        n_distro_cen_11A[i]=0;
        n_distro_cen_11B[i]=0;
        n_distro_cen_11C[i]=0;
        n_distro_cen_11W[i]=0;
        n_distro_cen_12A[i]=0;
        n_distro_cen_12B[i]=0;
        n_distro_cen_12K[i]=0;
        n_distro_cen_13A[i]=0;
        n_distro_cen_13B[i]=0;
        n_distro_cen_13K[i]=0;
        n_distro_cen_FCC[i]=0;
        n_distro_cen_HCP[i]=0;
        n_distro_cen_BCC_9[i]=0;
        n_distro_cen_BCC_15[i]=0;
    }
    
    for (i=0; i<9; i++) n_distro_shell_9B[i]=0;
    for (i=0; i<9; i++) n_distro_shell_9K[i]=0;
    for (i=0; i<10; i++) n_distro_shell_10B[i]=0;
    for (i=0; i<10; i++) n_distro_shell_10K[i]=0;
    for (i=0; i<10; i++) n_distro_shell_10W[i]=0;
    for (i=0; i<11; i++) n_distro_shell_11A[i]=0;
    for (i=0; i<11; i++) n_distro_shell_11B[i]=0;
    for (i=0; i<11; i++) n_distro_shell_11C[i]=0;
    for (i=0; i<11; i++) n_distro_shell_11W[i]=0;
    for (i=0; i<12; i++) n_distro_shell_12A[i]=0;
    for (i=0; i<12; i++) n_distro_shell_12B[i]=0;
    for (i=0; i<12; i++) n_distro_shell_12K[i]=0;
    for (i=0; i<13; i++) n_distro_shell_13A[i]=0;
    for (i=0; i<13; i++) n_distro_shell_13B[i]=0;
    for (i=0; i<13; i++) n_distro_shell_13K[i]=0;
    for (i=0; i<13; i++) n_distro_shell_FCC[i]=0;
    for (i=0; i<13; i++) n_distro_shell_HCP[i]=0;
    for (i=0; i<9; i++) n_distro_shell_BCC_9[i]=0;
    for (i=0; i<15; i++) n_distro_shell_BCC_15[i]=0;
    
    nAsp3=nAsp3a=nAsp3b=nAsp3c=0;
    nAsp4=nAsp4a=nAsp4b=nAsp4c=nA6A=0;
    nAsp5=nAsp5a=nAsp5b=nAsp5c=0;
    nA6Z=nA7K=0;
    nA8A=nA8B=nA8K=0;
    nA9A=nA9B=nA9K=0;
    nA10A=nA10B=nA10K=nA10W=0;
    nA11A=nA11B=nA11C=nA11E=nA11F=nA11W=0;
    nA12A=nA12B=nA12D=nA12E=nA12K=0;
    nA13A=nA13B=nA13K=0;
    nAFCC=nAHCP=nABCC_9=nABCC_15=0;

    nA_cen_9B=nA_cen_9K=0;
    nA_cen_10B=nA_cen_10K=nA_cen_10W=0;
    nA_cen_11A=nA_cen_11B=nA_cen_11C=nA_cen_11W=0;
    nA_cen_12A=nA_cen_12B=nA_cen_12K=0;
    nA_cen_13A=nA_cen_13B=nA_cen_13K=0;
    nA_cen_FCC=nA_cen_HCP=nA_cen_BCC_9=nA_cen_BCC_15=0;

    nA_shell_9B=nA_shell_9K=0;
    nA_shell_10B=nA_shell_10K=nA_shell_10W=0;
    nA_shell_11A=nA_shell_11B=nA_shell_11C=nA_shell_11W=0;
    nA_shell_12A=nA_shell_12B=nA_shell_12K=0;
    nA_shell_13A=nA_shell_13B=nA_shell_13K=0;
    nA_shell_FCC=nA_shell_HCP=nA_shell_BCC_9=nA_shell_BCC_15=0;

    nBsp3=nBsp3a=nBsp3b=nBsp3c=0;
    nBsp4=nBsp4a=nBsp4b=nBsp4c=nB6A=0;
    nBsp5=nBsp5a=nBsp5b=nBsp5c=0;
    nB6Z=nB7K=0;
    nB8A=nB8B=nB8K=0;
    nB9A=nB9B=nB9K=0;
    nB10A=nB10B=nB10K=nB10W=0;
    nB11A=nB11B=nB11C=nB11E=nB11F=nB11W=0;
    nB12A=nB12B=nB12D=nB12E=nB12K=0;
    nB13A=nB13B=nB13K=0;
    nBFCC=nBHCP=nBBCC_9=nBBCC_15=0;

    nB_cen_9B=nB_cen_9K=0;
    nB_cen_10B=nB_cen_10K=nB_cen_10W=0;
    nB_cen_11A=nB_cen_11B=nB_cen_11C=nB_cen_11W=0;
    nB_cen_12A=nB_cen_12B=nB_cen_12K=0;
    nB_cen_13A=nB_cen_13B=nB_cen_13K=0;
    nB_cen_FCC=nB_cen_HCP=nB_cen_BCC_9=nB_cen_BCC_15=0;

    nB_shell_9B=nB_shell_9K=0;
    nB_shell_10B=nB_shell_10K=nB_shell_10W=0;
    nB_shell_11A=nB_shell_11B=nB_shell_11C=nB_shell_11W=0;
    nB_shell_12A=nB_shell_12B=nB_shell_12K=0;
    nB_shell_13A=nB_shell_13B=nB_shell_13K=0;
    nB_shell_FCC=nB_shell_HCP=nB_shell_BCC_9=nB_shell_BCC_15=0;
}

void Setup_InitgsblVars(char *filename) { // Initialize ground state bond length deviation distribution arrays
    int i;
    char errMsg[1000];
    FILE *fin;
    
    printf("d%d reading bond length parameters from %s\n",rank,filename);
    fin=fopen(filename,"r");
    if (fin==NULL)  {
        sprintf(errMsg,"Setup_InitgsblVars() : Error opening file %s",filename);    // Always test file open
        Error_no_free(errMsg);
    }
    
    gsbl_sp3=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp3a=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp3b=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp3c=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_sp3 %lg gsbl_sp3a %lg gsbl_sp3b %lg gsbl_sp3c %lg\n",rank,gsbl_sp3,gsbl_sp3a,gsbl_sp3b,gsbl_sp3c);
    gsbl_sp4=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp4a=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp4b=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp4c=Setup_GetFirstDoubleFromLine(fin);
    gsbl_6A=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_sp4 %lg gsbl_sp4a %lg gsbl_sp4b %lg gsbl_sp4c %lg gsbl_6A %lg\n",rank,gsbl_sp4,gsbl_sp4a,gsbl_sp4b,gsbl_sp4c,gsbl_6A);
    gsbl_sp5=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp5a=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp5b=Setup_GetFirstDoubleFromLine(fin);
    gsbl_sp5c=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_sp5 %lg gsbl_sp5a %lg gsbl_sp5b %lg gsbl_sp5c %lg\n",rank,gsbl_sp5,gsbl_sp5a,gsbl_sp5b,gsbl_sp5c);
    gsbl_6Z=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_6Z %lg\n",rank,gsbl_6Z);
    gsbl_7K=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_7K %lg\n",rank,gsbl_7K);
    gsbl_8A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_8B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_8K=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_8A %lg gsbl_8B %lg gsbl_8K %lg\n",rank,gsbl_8A,gsbl_8B,gsbl_8K);
    gsbl_9A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_9B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_9K=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_9A %lg gsbl_9B %lg gsbl_9K %lg\n",rank,gsbl_9A,gsbl_9B,gsbl_9K);
    gsbl_10A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_10B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_10K=Setup_GetFirstDoubleFromLine(fin);
    gsbl_10W=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_10A %lg gsbl_10B %lg gsbl_10K %lg gsbl_10W %lg\n",rank,gsbl_10A,gsbl_10B,gsbl_10K,gsbl_10W);
    gsbl_11A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_11B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_11C=Setup_GetFirstDoubleFromLine(fin);
    gsbl_11E=Setup_GetFirstDoubleFromLine(fin);
    gsbl_11F=Setup_GetFirstDoubleFromLine(fin);
    gsbl_11W=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_11A %lg gsbl_11B %lg gsbl_11C %lg gsbl_11E %lg gsbl_11F %lg gsbl_11W %lg\n",rank,gsbl_11A,gsbl_11B,gsbl_11C,gsbl_11E,gsbl_11F,gsbl_11W);
    gsbl_12A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_12B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_12D=Setup_GetFirstDoubleFromLine(fin);
    gsbl_12E=Setup_GetFirstDoubleFromLine(fin);
    gsbl_12K=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_12A %lg gsbl_12B %lg gsbl_12D %lg gsbl_12E %lg gsbl_12K %lg\n",rank,gsbl_12A,gsbl_12B,gsbl_12D,gsbl_12E,gsbl_12K);
    gsbl_13A=Setup_GetFirstDoubleFromLine(fin);
    gsbl_13B=Setup_GetFirstDoubleFromLine(fin);
    gsbl_13K=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_13A %lg gsbl_13B %lg gsbl_13K %lg\n",rank,gsbl_13A,gsbl_13B,gsbl_13K);
    gsbl_FCC=Setup_GetFirstDoubleFromLine(fin);
    gsbl_HCP=Setup_GetFirstDoubleFromLine(fin);
    gsbl_BCC_9=Setup_GetFirstDoubleFromLine(fin);
    gsbl_BCC_15=Setup_GetFirstDoubleFromLine(fin);
    printf("d%d gsbl_FCC %lg gsbl_HCP %lg gsbl_BCC_9 %lg gsbl_BCC_15 %lg\n\n",rank,gsbl_FCC,gsbl_HCP,gsbl_BCC_9,gsbl_BCC_15);
    
    fclose(fin);
    
    bl_mom_sp3 = malloc(msp3*sizeof(double));   if (bl_mom_sp3==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp3[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_sp3a = malloc(msp3a*sizeof(double)); if (bl_mom_sp3a==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp3a[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp3b = malloc(msp3b*sizeof(double)); if (bl_mom_sp3b==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp3b[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp3c = malloc(msp3c*sizeof(double)); if (bl_mom_sp3c==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp3c[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp4 = malloc(msp4*sizeof(double));   if (bl_mom_sp4==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp4[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_sp4a = malloc(msp4a*sizeof(double)); if (bl_mom_sp4a==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp4a[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp4b = malloc(msp4b*sizeof(double)); if (bl_mom_sp4b==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp4b[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp4c = malloc(msp4c*sizeof(double)); if (bl_mom_sp4c==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp4c[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_6A = malloc(m6A*sizeof(double)); if (bl_mom_6A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_6A[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp5 = malloc(msp5*sizeof(double));   if (bl_mom_sp5==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp5[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_sp5a = malloc(msp5a*sizeof(double)); if (bl_mom_sp5a==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp5a[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp5b = malloc(msp5b*sizeof(double)); if (bl_mom_sp5b==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp5b[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_sp5c = malloc(msp5c*sizeof(double)); if (bl_mom_sp5c==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_sp5c[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    
    bl_mom_6Z = malloc(m6Z*sizeof(double)); if (bl_mom_6Z==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_6Z[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    
    bl_mom_7K = malloc(m7K*sizeof(double)); if (bl_mom_7K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_7K[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }

    bl_mom_8A = malloc(m8A*sizeof(double)); if (bl_mom_8A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_8A[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_8B = malloc(m8B*sizeof(double)); if (bl_mom_8B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_8B[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_8K = malloc(m8K*sizeof(double)); if (bl_mom_8K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_8K[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    
    bl_mom_9A = malloc(m9A*sizeof(double)); if (bl_mom_9A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_9A[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_9B = malloc(m9B*sizeof(double)); if (bl_mom_9B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_9B[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    bl_mom_9K = malloc(m9K*sizeof(double)); if (bl_mom_9K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_9K[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    
    bl_mom_10A = malloc(m10A*sizeof(double));   if (bl_mom_10A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_10A[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_10B = malloc(m10B*sizeof(double));   if (bl_mom_10B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_10B[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_10K = malloc(m10K*sizeof(double));   if (bl_mom_10K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_10K[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_10W = malloc(m10W*sizeof(double));   if (bl_mom_10W==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_10W[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    
    bl_mom_11A = malloc(m11A*sizeof(double));   if (bl_mom_11A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11A[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_11B = malloc(m11B*sizeof(double));   if (bl_mom_11B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11B[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_11C = malloc(m11C*sizeof(double));   if (bl_mom_11C==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11C[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_11E = malloc(m11E*sizeof(double));   if (bl_mom_11E==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11E[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_11F = malloc(m11F*sizeof(double));   if (bl_mom_11F==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11F[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_11W = malloc(m11W*sizeof(double));   if (bl_mom_11W==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_11W[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    
    bl_mom_12A = malloc(m12A*sizeof(double));   if (bl_mom_12A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_12A[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_12B = malloc(m12B*sizeof(double));   if (bl_mom_12B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_12B[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_12D = malloc(m12D*sizeof(double));   if (bl_mom_12D==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_12D[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_12E = malloc(m12E*sizeof(double));   if (bl_mom_12E==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_12E[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_12K = malloc(m12K*sizeof(double));   if (bl_mom_12K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_12K[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    
    bl_mom_13A = malloc(m13A*sizeof(double));   if (bl_mom_13A==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_13A[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_13B = malloc(m13B*sizeof(double));   if (bl_mom_13B==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_13B[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_13K = malloc(m13K*sizeof(double));   if (bl_mom_13K==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_13K[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    
    bl_mom_FCC = malloc(mFCC*sizeof(double));   if (bl_mom_FCC==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_FCC[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_HCP = malloc(mHCP*sizeof(double));   if (bl_mom_HCP==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_HCP[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_BCC_9 = malloc(mBCC_9*sizeof(double));   if (bl_mom_BCC_9==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_BCC_9[] malloc out of memory\bl_mom_");    Error_no_free(errMsg); }
    bl_mom_BCC_15 = malloc(mBCC_15*sizeof(double)); if (bl_mom_BCC_15==NULL) { sprintf(errMsg,"Setup_Ibl_mom_itgsblVars(): bl_mom_BCC_15[] malloc out of memory\bl_mom_");  Error_no_free(errMsg); }
    
    for (i=0; i<msp3; ++i) bl_mom_sp3[i]=0.0;
    for (i=0; i<msp3a; ++i) bl_mom_sp3a[i]=0.0;
    for (i=0; i<msp3b; ++i) bl_mom_sp3b[i]=0.0;
    for (i=0; i<msp3c; ++i) bl_mom_sp3c[i]=0.0;
    for (i=0; i<msp4; ++i) bl_mom_sp4[i]=0.0;
    for (i=0; i<msp4a; ++i) bl_mom_sp4a[i]=0.0;
    for (i=0; i<msp4b; ++i) bl_mom_sp4b[i]=0.0;
    for (i=0; i<msp4c; ++i) bl_mom_sp4c[i]=0.0;
    for (i=0; i<m6A; ++i) bl_mom_6A[i]=0.0;
    for (i=0; i<msp5; ++i) bl_mom_sp5[i]=0.0;
    for (i=0; i<msp5a; ++i) bl_mom_sp5a[i]=0.0;
    for (i=0; i<msp5b; ++i) bl_mom_sp5b[i]=0.0;
    for (i=0; i<msp5c; ++i) bl_mom_sp5c[i]=0.0;
    for (i=0; i<m6Z; ++i) bl_mom_6Z[i]=0.0;
    for (i=0; i<m7K; ++i) bl_mom_7K[i]=0.0;
    for (i=0; i<m8A; ++i) bl_mom_8A[i]=0.0;
    for (i=0; i<m8B; ++i) bl_mom_8B[i]=0.0;
    for (i=0; i<m8K; ++i) bl_mom_8K[i]=0.0;
    for (i=0; i<m9A; ++i) bl_mom_9A[i]=0.0;
    for (i=0; i<m9B; ++i) bl_mom_9B[i]=0.0;
    for (i=0; i<m9K; ++i) bl_mom_9K[i]=0.0;
    for (i=0; i<m10A; ++i) bl_mom_10A[i]=0.0;
    for (i=0; i<m10B; ++i) bl_mom_10B[i]=0.0;
    for (i=0; i<m10K; ++i) bl_mom_10K[i]=0.0;
    for (i=0; i<m10W; ++i) bl_mom_10W[i]=0.0;
    for (i=0; i<m11A; ++i) bl_mom_11A[i]=0.0;
    for (i=0; i<m11B; ++i) bl_mom_11B[i]=0.0;
    for (i=0; i<m11C; ++i) bl_mom_11C[i]=0.0;
    for (i=0; i<m11E; ++i) bl_mom_11E[i]=0.0;
    for (i=0; i<m11F; ++i) bl_mom_11F[i]=0.0;
    for (i=0; i<m11W; ++i) bl_mom_11W[i]=0.0;
    for (i=0; i<m12A; ++i) bl_mom_12A[i]=0.0;
    for (i=0; i<m12B; ++i) bl_mom_12B[i]=0.0;
    for (i=0; i<m12D; ++i) bl_mom_12D[i]=0.0;
    for (i=0; i<m12E; ++i) bl_mom_12E[i]=0.0;
    for (i=0; i<m12K; ++i) bl_mom_12K[i]=0.0;
    for (i=0; i<m13A; ++i) bl_mom_13A[i]=0.0;
    for (i=0; i<m13B; ++i) bl_mom_13B[i]=0.0;
    for (i=0; i<m13K; ++i) bl_mom_13K[i]=0.0;
    for (i=0; i<mFCC; ++i) bl_mom_FCC[i]=0.0;
    for (i=0; i<mHCP; ++i) bl_mom_HCP[i]=0.0;
    for (i=0; i<mBCC_9; ++i) bl_mom_BCC_9[i]=0.0;
    for (i=0; i<mBCC_15; ++i) bl_mom_BCC_15[i]=0.0;

    mean_bl_mom_sp3=mean_bl_mom_sp3a=mean_bl_mom_sp3b=mean_bl_mom_sp3c=0.0; // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
    mean_bl_mom_sp4=mean_bl_mom_sp4a=mean_bl_mom_sp4b=mean_bl_mom_sp4c=0.0; // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
    mean_bl_mom_sp5=mean_bl_mom_sp5a=mean_bl_mom_sp5b=mean_bl_mom_sp5c=0.0; // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
    mean_bl_mom_6A=mean_bl_mom_6Z=mean_bl_mom_7K=0.0;   // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_8A=mean_bl_mom_8B=mean_bl_mom_8K=0.0;   // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_9A=mean_bl_mom_9B=mean_bl_mom_9K=0.0;   // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_10A=mean_bl_mom_10B=mean_bl_mom_10K=mean_bl_mom_10W=0.0;    // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_11A=mean_bl_mom_11B=mean_bl_mom_11C=mean_bl_mom_11E=mean_bl_mom_11F=mean_bl_mom_11W=0.0;    // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_12A=mean_bl_mom_12B=mean_bl_mom_12D=mean_bl_mom_12E=mean_bl_mom_12K=0.0;    // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_13A=mean_bl_mom_13B=mean_bl_mom_13K=0.0;    // mean_bl_mom_ax size of mean_bl_mom_** arrays in dimean_bl_mom_ension i
    mean_bl_mom_FCC=mean_bl_mom_HCP=mean_bl_mom_BCC_9=mean_bl_mom_BCC_15=0.0;   // mean_bl_mom_ax size of **sp** arrays in dimean_bl_mom_ension i
}
void Setup_InitDynamicVars(char *filename) { // Initialize dynamic variables
    int i,j;
    char errMsg[1000];
    FILE *fin;
    
    printf("d%d reading dynamics memory parameters from %s\n",rank,filename);
    fin=fopen(filename,"r");
    if (fin==NULL)  {
        sprintf(errMsg,"Setup_InitDynamicVars() : Error opening file %s",filename); // Always test file open
        Error_no_free(errMsg);
    }
    
    dosp3=dosp3a=dosp3b=dosp3c=0;
    dosp4=dosp4a=dosp4b=dosp4c=0;
    dosp5=dosp5a=dosp5b=dosp5c=0;
    do6Z=do7K=do8A=do8B=do8K=do9A=do9B=do9K=do10A=do10B=do10K=do10W=0;
    do11A=do11B=do11C=do11E=do11F=do11W=do12A=do12B=do12D=do12E=do12K=0;
    do13A=do13B=do13K=doFCC=doHCP=doBCC9=doBCC15=0;
    
    dyn_MaxLives=dyn_TotClus=dyn_TotEvents=dyn_Tot_a_Events=0;
    
    dyn_msp3=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp3>0) { dyn_msp3=initNoDynamicClusters; dosp3=1; dosp3a=1; dosp3b=1; dosp3c=1; }
    dyn_msp3a=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp3a>0) { dyn_msp3a=initNoDynamicClusters; dosp3=1; dosp3a=1; dosp3b=1; dosp3c=1; }
    dyn_msp3b=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp3b>0) { dyn_msp3b=initNoDynamicClusters; dosp3=1; dosp3a=1; dosp3b=1; dosp3c=1; }
    dyn_msp3c=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp3c>0) { dyn_msp3c=initNoDynamicClusters; dosp3=1; dosp3a=1; dosp3b=1; dosp3c=1; }
    printf("d%d dyn_msp3 %d dyn_msp3a %d dyn_msp3b %d dyn_msp3c %d\n",rank,dyn_msp3,dyn_msp3a,dyn_msp3b,dyn_msp3c);
    dyn_msp4=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp4>0) { dyn_msp4=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4a=1; dosp4b=1; dosp4c=1; }
    dyn_msp4a=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp4a>0) { dyn_msp4a=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4a=1; dosp4b=1; dosp4c=1; } 
    dyn_msp4b=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp4b>0) { dyn_msp4b=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4a=1; dosp4b=1; dosp4c=1; }
    dyn_m6A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m6A>0) { dyn_m6A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4a=1; dosp4b=1; dosp4c=1; }
    printf("d%d dyn_msp4 %d dyn_msp4a %d dyn_msp4b %d dyn_m6A %d\n",rank,dyn_msp4,dyn_msp4a,dyn_msp4b,dyn_m6A);
    dyn_msp5=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp5>0) { dyn_msp5=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5a=1; dosp5b=1; dosp5c=1; }
    dyn_msp5a=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp5a>0) { dyn_msp5a=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5a=1; dosp5b=1; dosp5c=1; }
    dyn_msp5b=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp5b>0) { dyn_msp5b=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5a=1; dosp5b=1; dosp5c=1; }
    dyn_msp5c=Setup_GetFirstIntFromLine(fin);
    if (dyn_msp5c>0) { dyn_msp5c=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5a=1; dosp5b=1; dosp5c=1; }
    printf("d%d dyn_msp5 %d dyn_msp5a %d dyn_msp5b %d dyn_msp5c %d\n",rank,dyn_msp5,dyn_msp5a,dyn_msp5b,dyn_msp5c);
    dyn_m6Z=Setup_GetFirstIntFromLine(fin);
    if (dyn_m6Z>0) { dyn_m6Z=initNoDynamicClusters; dosp3=1; dosp3c=1;  do6Z=1; }
    printf("d%d dyn_m6Z %d\n",rank,dyn_m6Z);
    dyn_m7K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m7K>0) { dyn_m7K=initNoDynamicClusters; dosp3=1; dosp3c=1;  do7K=1; }
    printf("d%d dyn_m7K %d\n",rank,dyn_m7K);
    dyn_m8A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m8A>0) { dyn_m8A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5b=1;  dosp5c=1; do8A=1; }
    dyn_m8B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m8B>0) { dyn_m8B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do8B=1; }
    dyn_m8K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m8K>0) { dyn_m8K=initNoDynamicClusters; dosp3=1; dosp3c=1;  do8K=1; }
    printf("d%d dyn_m8A %d dyn_m8B %d dyn_m8K %d\n",rank,dyn_m8A,dyn_m8B,dyn_m8K);
    dyn_m9A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m9A>0) { dyn_m9A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4b=1; do9A=1; }
    dyn_m9B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m9B>0) { dyn_m9B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; }
    dyn_m9K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m9K>0) { dyn_m9K=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4c=1; do9K=1; }
    printf("d%d dyn_m9A %d dyn_m9B %d dyn_m9K %d\n",rank,dyn_m9A,dyn_m9B,dyn_m9K);
    dyn_m10A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m10A>0) { dyn_m10A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4b=1; do10A=1; }
    dyn_m10B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m10B>0) { dyn_m10B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; do10B=1; }
    dyn_m10K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m10K>0) { dyn_m10K=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9K=1; do10K=1; }
    dyn_m10W=Setup_GetFirstIntFromLine(fin);
    if (dyn_m10W>0) { dyn_m10W=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5b=1; do10W=1; }
    printf("d%d dyn_m10A %d dyn_m10B %d dyn_m10K %d dyn_m10W %d\n",rank,dyn_m10A,dyn_m10B,dyn_m10K,dyn_m10W);
    dyn_m11A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11A>0) { dyn_m11A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4c=1; do11A=1; }
    dyn_m11B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11B>0) { dyn_m11B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; do11B=1; }
    dyn_m11C=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11C>0) { dyn_m11C=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do11C=1; }
    dyn_m11E=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11E>0) { dyn_m11E=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; do11E=1; }
    dyn_m11F=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11F>0) { dyn_m11F=initNoDynamicClusters; dosp3=1; dosp3c=1; dosp4=1; dosp4c=1; do11F=1; }
    dyn_m11W=Setup_GetFirstIntFromLine(fin);
    if (dyn_m11W>0) { dyn_m11W=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; do10B=1; do11W=1; }
    printf("d%d dyn_m11A %d dyn_m11B %d dyn_m11C %d dyn_m11E %d dyn_m11F %d dyn_m11W %d\n",rank,dyn_m11A,dyn_m11B,dyn_m11C,dyn_m11E,dyn_m11F,dyn_m11W);
    dyn_m12A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m12A>0) { dyn_m12A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do11C=1; do12A=1; }
    dyn_m12B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m12B>0) { dyn_m12B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do12B=1; }
    dyn_m12D=Setup_GetFirstIntFromLine(fin);
    if (dyn_m12D>0) { dyn_m12D=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do9B=1; do11E=1; do12D=1; }
    dyn_m12E=Setup_GetFirstIntFromLine(fin);
    if (dyn_m12E>0) { dyn_m12E=initNoDynamicClusters; dosp3=1; dosp3c=1; dosp4=1; dosp4c=1; do11F=1; do12E=1; }
    dyn_m12K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m12K>0) { dyn_m12K=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4c=1; do11A=1; do12K=1; }
    printf("d%d dyn_m12A %d dyn_m12B %d dyn_m12D %d dyn_m12E %d dyn_m12K %d\n",rank,dyn_m12A,dyn_m12B,dyn_m12D,dyn_m12E,dyn_m12K);
    dyn_m13A=Setup_GetFirstIntFromLine(fin);
    if (dyn_m13A>0) { dyn_m13A=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do12B=1;  do13A=1;}
    dyn_m13B=Setup_GetFirstIntFromLine(fin);
    if (dyn_m13B>0) { dyn_m13B=initNoDynamicClusters; dosp3=1; dosp4=1; dosp5=1; dosp5c=1; do13B=1; }
    dyn_m13K=Setup_GetFirstIntFromLine(fin);
    if (dyn_m13K>0) { dyn_m13K=initNoDynamicClusters; dosp3=1; dosp3c=1; dosp4=1; dosp4c=1; do11F=1; do13K=1; }
    printf("d%d dyn_m13A %d dyn_m13B %d dyn_m13K %d\n",rank,dyn_m13A,dyn_m13B,dyn_m13K);
    dyn_mFCC=Setup_GetFirstIntFromLine(fin);
    if (dyn_mFCC>0) { dyn_mFCC=initNoDynamicClusters; dosp3=1; dosp3b=1; dosp3c=1; doFCC=1; }
    dyn_mHCP=Setup_GetFirstIntFromLine(fin);
    if (dyn_mHCP>0) { dyn_mHCP=initNoDynamicClusters; dosp3=1; dosp3c=1; doHCP=1; }
    dyn_mBCC_9=Setup_GetFirstIntFromLine(fin);
    if (dyn_mBCC_9>0) { dyn_mBCC_9=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4b=1; dosp4c=1; doBCC9=1; }
    dyn_mBCC_15=Setup_GetFirstIntFromLine(fin);
    if (dyn_mBCC_15>0) { dyn_mBCC_15=initNoDynamicClusters; dosp3=1; dosp4=1; dosp4c=1; doBCC15=1; }
    printf("d%d dyn_mFCC %d dyn_mHCP %d dyn_mBCC_9 %d dyn_mBCC_15 %d\n\n",rank,dyn_mFCC,dyn_mHCP,dyn_mBCC_9,dyn_mBCC_15);
    
    printf("d%d dosp3 %d dosp3a %d dosp3b %d dosp3c %d\n",rank,dosp3,dosp3a,dosp3b,dosp3c);
    printf("d%d dosp4 %d dosp4a %d dosp4b %d dosp4c %d\n",rank,dosp4,dosp4a,dosp4b,dosp4c);
    printf("d%d dosp5 %d dosp5a %d dosp5b %d dosp5c %d\n",rank,dosp5,dosp5a,dosp5b,dosp5c);
    printf("d%d do6Z %d do7K %d do8A %d do8B %d do8K %d do9A %d do9B %d do9K %d do10A %d do10B %d do10K %d do10W %d\n",rank,do6Z,do7K,do8A,do8B,do8K,do9A,do9B,do9K,do10A,do10B,do10K,do10W);
    printf("d%d do11A %d do11B %d do11C %d do11E %d do11F %d do11W %d do12A %d do12B %d do12D %d do12E %d do12K %d\n",rank,do11A,do11B,do11C,do11E,do11F,do11W,do12A,do12B,do12D,do12E,do12K);
    printf("d%d do13A %d do13B %d do13K %d doFCC %d doHCP %d doBCC9 %d doBCC15 %d\n\n",rank,do13A,do13B,do13K,doFCC,doHCP,doBCC9,doBCC15);
    
    fclose(fin);
    
    dyn_nsp3=dyn_nsp3a=dyn_nsp3b=dyn_nsp3c=0;
    dyn_nsp4=dyn_nsp4a=dyn_nsp4b=dyn_n6A=0;
    dyn_nsp5=dyn_nsp5a=dyn_nsp5b=dyn_nsp5c=0;
    
    dyn_n6Z=dyn_n7K=0;
    dyn_n8A=dyn_n8B=dyn_n8K=0;
    dyn_n9A=dyn_n9B=dyn_n9K=0;
    dyn_n10A=dyn_n10B=dyn_n10K=dyn_n10W=0;
    dyn_n11A=dyn_n11B=dyn_n11C=dyn_n11E=dyn_n11F=dyn_n11W=0;
    dyn_n12A=dyn_n12B=dyn_n12D=dyn_n12E=dyn_n12K=0;
    dyn_n13A=dyn_n13B=dyn_n13K=0;
    dyn_nFCC=dyn_nHCP=dyn_nBCC_9=dyn_nBCC_15=0;
    
    dyn_esp3=dyn_esp3a=dyn_esp3b=dyn_esp3c=0;
    dyn_esp4=dyn_esp4a=dyn_esp4b=dyn_e6A=0;
    dyn_esp5=dyn_esp5a=dyn_esp5b=dyn_esp5c=0;
    
    dyn_e6Z=dyn_e7K=0;
    dyn_e8A=dyn_e8B=dyn_e8K=0;
    dyn_e9A=dyn_e9B=dyn_e9K=0;
    dyn_e10A=dyn_e10B=dyn_e10K=dyn_e10W=0;
    dyn_e11A=dyn_e11B=dyn_e11C=dyn_e11E=dyn_e11F=dyn_e11W=0;
    dyn_e12A=dyn_e12B=dyn_e12D=dyn_e12E=dyn_e12K=0;
    dyn_e13A=dyn_e13B=dyn_e13K=0;
    dyn_eFCC=dyn_eHCP=dyn_eBCC_9=dyn_eBCC_15=0;
    
    // dynamic clusters
    if (dyn_msp3!=-1) {
        dyn_sp3 = malloc(dyn_msp3*sizeof(int *));   if (dyn_sp3==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3; ++j) { dyn_sp3[j] = malloc(3*sizeof(int));    if (dyn_sp3[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lsp3 = malloc(dyn_msp3*sizeof(int *));  if (dyn_lsp3==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3; ++j) { dyn_lsp3[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lsp3[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp3; ++i) {
            for (j=0;j<3;++j) dyn_sp3[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp3[i][j]=-1;
        }
    }
    if (dyn_msp3a!=-1) {
        dyn_sp3a = malloc(dyn_msp3a*sizeof(int *)); if (dyn_sp3a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3a; ++j) { dyn_sp3a[j] = malloc(3*sizeof(int));  if (dyn_sp3a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3a[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp3a = malloc(dyn_msp3a*sizeof(int *));    if (dyn_lsp3a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3a[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3a; ++j) { dyn_lsp3a[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp3a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3a[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp3a; ++i) {
            for (j=0;j<3;++j) dyn_sp3a[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp3a[i][j]=-1;
        }
    }
    if (dyn_msp3b!=-1) {
        dyn_sp3b = malloc(dyn_msp3b*sizeof(int *)); if (dyn_sp3b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3b; ++j) { dyn_sp3b[j] = malloc(4*sizeof(int));  if (dyn_sp3b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3b[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp3b = malloc(dyn_msp3b*sizeof(int *));    if (dyn_lsp3b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3b[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3b; ++j) { dyn_lsp3b[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp3b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3b[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp3b = malloc(msp3b*sizeof(int));    if (dyn_up_sp3b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp3b[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_msp3b; ++i) {
            for (j=0;j<4;++j) dyn_sp3b[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp3b[i][j]=-1;
        }
        for (i=0; i<msp3b; ++i) {
            if (doSubClusts==1) dyn_up_sp3b[i]=-1;
        }
    } 
    if (dyn_msp3c!=-1) {
        dyn_sp3c = malloc(dyn_msp3c*sizeof(int *)); if (dyn_sp3c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3c[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3c; ++j) { dyn_sp3c[j] = malloc(5*sizeof(int));  if (dyn_sp3c[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp3c[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp3c = malloc(dyn_msp3c*sizeof(int *));    if (dyn_lsp3c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3c[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp3c; ++j) { dyn_lsp3c[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp3c[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp3c[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp3c = malloc(msp3c*sizeof(int));    if (dyn_up_sp3c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp3c[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_msp3c; ++i) {
            for (j=0;j<5;++j) dyn_sp3c[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp3c[i][j]=-1;
        }
        for (i=0; i<msp3c; ++i) {
            if (doSubClusts==1) dyn_up_sp3c[i]=-1;
        }   
    }
    if (dyn_msp4!=-1) {
        dyn_sp4 = malloc(dyn_msp4*sizeof(int *));   if (dyn_sp4==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4; ++j) { dyn_sp4[j] = malloc(4*sizeof(int));    if (dyn_sp4[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lsp4 = malloc(dyn_msp4*sizeof(int *));  if (dyn_lsp4==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4; ++j) { dyn_lsp4[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lsp4[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp4; ++i) {
            for (j=0;j<4;++j) dyn_sp4[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp4[i][j]=-1;
        }
    }
    if (dyn_msp4a!=-1) {
        dyn_sp4a = malloc(dyn_msp4a*sizeof(int *)); if (dyn_sp4a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4a; ++j) { dyn_sp4a[j] = malloc(4*sizeof(int));  if (dyn_sp4a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4a[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp4a = malloc(dyn_msp4a*sizeof(int *));    if (dyn_lsp4a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4a[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4a; ++j) { dyn_lsp4a[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp4a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4a[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp4a; ++i) {
            for (j=0;j<4;++j) dyn_sp4a[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp4a[i][j]=-1;
        }
    }
    if (dyn_msp4b!=-1) {
        dyn_sp4b = malloc(dyn_msp4b*sizeof(int *)); if (dyn_sp4b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4b; ++j) { dyn_sp4b[j] = malloc(5*sizeof(int));  if (dyn_sp4b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp4b[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp4b = malloc(dyn_msp4b*sizeof(int *));    if (dyn_lsp4b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4b[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp4b; ++j) { dyn_lsp4b[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp4b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp4b[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp4b = malloc(msp4b*sizeof(int));    if (dyn_up_sp4b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp4b[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_msp4b; ++i) {
            for (j=0;j<5;++j) dyn_sp4b[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp4b[i][j]=-1;
        }
        for (i=0; i<msp4b; ++i) {
            if (doSubClusts==1) dyn_up_sp4b[i]=-1;
        }
    }
    if (dyn_m6A!=-1) {
        dyn_hc6A = malloc(dyn_m6A*sizeof(int *));   if (dyn_hc6A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc6A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m6A; ++j) { dyn_hc6A[j] = malloc(6*sizeof(int));    if (dyn_hc6A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc6A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l6A = malloc(dyn_m6A*sizeof(int *));    if (dyn_l6A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l6A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m6A; ++j) { dyn_l6A[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l6A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l6A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp4c = malloc(msp4c*sizeof(int));    if (dyn_up_sp4c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp4c[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m6A; ++i) {
            for (j=0;j<6;++j) dyn_hc6A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l6A[i][j]=-1;
        }
        for (i=0; i<msp4c; ++i) {
            if (doSubClusts==1) dyn_up_sp4c[i]=-1;
        }
    }
    if (dyn_msp5!=-1) {
        dyn_sp5 = malloc(dyn_msp5*sizeof(int *));   if (dyn_sp5==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5; ++j) { dyn_sp5[j] = malloc(5*sizeof(int));    if (dyn_sp5[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lsp5 = malloc(dyn_msp5*sizeof(int *));  if (dyn_lsp5==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5; ++j) { dyn_lsp5[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lsp5[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp5; ++i) { 
            for (j=0;j<5;++j) dyn_sp5[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp5[i][j]=-1;
        }
    }
    if (dyn_msp5a!=-1) {
        dyn_sp5a = malloc(dyn_msp5a*sizeof(int *)); if (dyn_sp5a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5a[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5a; ++j) { dyn_sp5a[j] = malloc(5*sizeof(int));  if (dyn_sp5a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5a[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp5a = malloc(dyn_msp5a*sizeof(int *));    if (dyn_lsp5a==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5a[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5a; ++j) { dyn_lsp5a[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp5a[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5a[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        for (i=0; i<dyn_msp5a; ++i) { 
            for (j=0;j<5;++j) dyn_sp5a[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp5a[i][j]=-1;
        }
    }
    if (dyn_msp5b!=-1) {
        dyn_sp5b = malloc(dyn_msp5b*sizeof(int *)); if (dyn_sp5b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5b[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5b; ++j) { dyn_sp5b[j] = malloc(6*sizeof(int));  if (dyn_sp5b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5b[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp5b = malloc(dyn_msp5b*sizeof(int *));    if (dyn_lsp5b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5b[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5b; ++j) { dyn_lsp5b[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp5b[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5b[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp5b = malloc(msp5b*sizeof(int));    if (dyn_up_sp5b==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp5b[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_msp5b; ++i) {
            for (j=0;j<6;++j) dyn_sp5b[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp5b[i][j]=-1;
        }
        for (i=0; i<msp5b; ++i) {
            if (doSubClusts==1) dyn_up_sp5b[i]=-1;
        }
    }
    if (dyn_msp5c!=-1) {
        dyn_sp5c = malloc(dyn_msp5c*sizeof(int *)); if (dyn_sp5c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5c[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5c; ++j) { dyn_sp5c[j] = malloc(7*sizeof(int));  if (dyn_sp5c[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sp5c[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lsp5c = malloc(dyn_msp5c*sizeof(int *));    if (dyn_lsp5c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5c[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_msp5c; ++j) { dyn_lsp5c[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lsp5c[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lsp5c[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1) {
            dyn_up_sp5c = malloc(msp5c*sizeof(int));    if (dyn_up_sp5c==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_sp5c[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_msp5c; ++i) {
            for (j=0;j<7;++j) dyn_sp5c[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lsp5c[i][j]=-1;
        }
        for (i=0; i<msp5c; ++i) {
            if (doSubClusts==1) dyn_up_sp5c[i]=-1;
        }   
    }
    if (dyn_m6Z!=-1) {
        dyn_hc6Z = malloc(dyn_m6Z*sizeof(int *));   if (dyn_hc6Z==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc6Z[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m6Z; ++j) { dyn_hc6Z[j] = malloc(6*sizeof(int));    if (dyn_hc6Z[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc6Z[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l6Z = malloc(dyn_m6Z*sizeof(int *));    if (dyn_l6Z==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l6Z[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m6Z; ++j) { dyn_l6Z[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l6Z[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l6Z[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            dyn_sub_6Z = malloc(dyn_m6Z*sizeof(int *)); if (dyn_sub_6Z==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_6Z[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m6Z; ++j) { dyn_sub_6Z[j] = malloc(2*sizeof(int));  if (dyn_sub_6Z[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_6Z[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m6Z; ++i) {
            for (j=0;j<6;++j) dyn_hc6Z[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l6Z[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1) for (j=0;j<2;++j) dyn_sub_6Z[i][j]=-1;
        }
    }
    if (dyn_m7K!=-1) {
        dyn_hc7K = malloc(dyn_m7K*sizeof(int *));   if (dyn_hc7K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc7K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m7K; ++j) { dyn_hc7K[j] = malloc(7*sizeof(int));    if (dyn_hc7K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc7K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l7K = malloc(dyn_m7K*sizeof(int *));    if (dyn_l7K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l7K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m7K; ++j) { dyn_l7K[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l7K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l7K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            dyn_sub_7K = malloc(dyn_m7K*sizeof(int *)); if (dyn_sub_7K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_7K[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m7K; ++j) { dyn_sub_7K[j] = malloc(2*sizeof(int));  if (dyn_sub_7K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_7K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m7K; ++i) {
            for (j=0;j<7;++j) dyn_hc7K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l7K[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1) for (j=0;j<2;++j) dyn_sub_7K[i][j]=-1;
        }
    }
    if (dyn_m8A!=-1) {
        dyn_hc8A = malloc(dyn_m8A*sizeof(int *));   if (dyn_hc8A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m8A; ++j) { dyn_hc8A[j] = malloc(8*sizeof(int));    if (dyn_hc8A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l8A = malloc(dyn_m8A*sizeof(int *));    if (dyn_l8A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m8A; ++j) { dyn_l8A[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l8A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
            dyn_sub_8A = malloc(dyn_m8A*sizeof(int *)); if (dyn_sub_8A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8A[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m8A; ++j) { dyn_sub_8A[j] = malloc(12*sizeof(int)); if (dyn_sub_8A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m8A; ++i){
            for (j=0;j<8;++j) dyn_hc8A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l8A[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) for (j=0;j<12;++j) dyn_sub_8A[i][j]=-1;
        }
    }
    if (dyn_m8B!=-1) {
        dyn_hc8B = malloc(dyn_m8B*sizeof(int *));   if (dyn_hc8B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m8B; ++j) { dyn_hc8B[j] = malloc(8*sizeof(int));    if (dyn_hc8B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l8B = malloc(dyn_m8B*sizeof(int *));    if (dyn_l8B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m8B; ++j) { dyn_l8B[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l8B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_8B = malloc(dyn_m8B*sizeof(int *)); if (dyn_sub_8B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8B[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m8B; ++j) { dyn_sub_8B[j] = malloc(1*sizeof(int));  if (dyn_sub_8B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m8B; ++i) {
            for (j=0;j<8;++j) dyn_hc8B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l8B[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<1;++j) dyn_sub_8B[i][j]=-1;
        }
    }
    if (dyn_m8K!=-1) {
        dyn_hc8K = malloc(dyn_m8K*sizeof(int *));   if (dyn_hc8K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m8K; ++j) { dyn_hc8K[j] = malloc(8*sizeof(int));    if (dyn_hc8K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc8K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l8K = malloc(dyn_m8K*sizeof(int *));    if (dyn_l8K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m8K; ++j) { dyn_l8K[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l8K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l8K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            dyn_sub_8K = malloc(dyn_m8K*sizeof(int *)); if (dyn_sub_8K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8K[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m8K; ++j) { dyn_sub_8K[j] = malloc(3*sizeof(int));  if (dyn_sub_8K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_8K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m8K; ++i) {
            for (j=0;j<8;++j) dyn_hc8K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l8K[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1) for (j=0;j<3;++j) dyn_sub_8K[i][j]=-1;
        }
    }
    if (dyn_m9A!=-1) {
        dyn_hc9A = malloc(dyn_m9A*sizeof(int *));   if (dyn_hc9A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m9A; ++j) { dyn_hc9A[j] = malloc(9*sizeof(int));    if (dyn_hc9A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l9A = malloc(dyn_m9A*sizeof(int *));    if (dyn_l9A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m9A; ++j) { dyn_l9A[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l9A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp4b!=-1) {
            dyn_sub_9A = malloc(dyn_m9A*sizeof(int *)); if (dyn_sub_9A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9A[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m9A; ++j) { dyn_sub_9A[j] = malloc(3*sizeof(int));  if (dyn_sub_9A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m9A; ++i){
            for (j=0;j<9;++j) dyn_hc9A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l9A[i][j]=-1;
            if (doSubClusts==1 && dyn_msp4b!=-1) for (j=0;j<3;++j) dyn_sub_9A[i][j]=-1;
        }
    }
    if (dyn_m9B!=-1) {
        dyn_hc9B = malloc(dyn_m9B*sizeof(int *));   if (dyn_hc9B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m9B; ++j) { dyn_hc9B[j] = malloc(9*sizeof(int));    if (dyn_hc9B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l9B = malloc(dyn_m9B*sizeof(int *));    if (dyn_l9B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m9B; ++j) { dyn_l9B[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l9B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_9B = malloc(dyn_m9B*sizeof(int *)); if (dyn_sub_9B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9B[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m9B; ++j) { dyn_sub_9B[j] = malloc(2*sizeof(int));  if (dyn_sub_9B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_9B = malloc(m9B*sizeof(int));    if (dyn_up_9B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_9B[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m9B; ++i) {
            for (j=0;j<9;++j) dyn_hc9B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l9B[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<2;++j) dyn_sub_9B[i][j]=-1;
        }
        for (i=0; i<m9B; ++i) {
            if (doSubClusts==1) dyn_up_9B[i]=-1;
        }
    }
    if (dyn_m9K!=-1) {
        dyn_hc9K = malloc(dyn_m9K*sizeof(int *));   if (dyn_hc9K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m9K; ++j) { dyn_hc9K[j] = malloc(9*sizeof(int));    if (dyn_hc9K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc9K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_l9K = malloc(dyn_m9K*sizeof(int *));    if (dyn_l9K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m9K; ++j) { dyn_l9K[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_l9K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l9K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m6A!=-1) {
            dyn_sub_9K = malloc(dyn_m9K*sizeof(int *)); if (dyn_sub_9K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9K[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_m9K; ++j) { dyn_sub_9K[j] = malloc(2*sizeof(int));  if (dyn_sub_9K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_9K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_9K = malloc(m9K*sizeof(int));    if (dyn_up_9K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_9K[] malloc out of memory\n");   Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m9K; ++i) {
            for (j=0;j<9;++j) dyn_hc9K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l9K[i][j]=-1;
            if (doSubClusts==1 && dyn_m6A!=-1) for (j=0;j<2;++j) dyn_sub_9K[i][j]=-1;
        }
        for (i=0; i<m9K; ++i) {
            if (doSubClusts==1) dyn_up_9K[i]=-1;
        }
    }
    if (dyn_m10A!=-1) {
        dyn_hc10A = malloc(dyn_m10A*sizeof(int *)); if (dyn_hc10A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m10A; ++j) { dyn_hc10A[j] = malloc(10*sizeof(int)); if (dyn_hc10A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l10A = malloc(dyn_m10A*sizeof(int *));  if (dyn_l10A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m10A; ++j) { dyn_l10A[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l10A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp4b!=-1) {
            dyn_sub_10A = malloc(dyn_m10A*sizeof(int *));   if (dyn_sub_10A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10A[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m10A; ++j) { dyn_sub_10A[j] = malloc(2*sizeof(int));    if (dyn_sub_10A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m10A; ++i) {
            for (j=0;j<10;++j) dyn_hc10A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l10A[i][j]=-1;
            if (doSubClusts==1 && dyn_msp4b!=-1) for (j=0;j<2;++j) dyn_sub_10A[i][j]=-1;
        }
    }
    if (dyn_m10B!=-1) {
        dyn_hc10B = malloc(dyn_m10B*sizeof(int *)); if (dyn_hc10B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m10B; ++j) { dyn_hc10B[j] = malloc(10*sizeof(int)); if (dyn_hc10B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l10B = malloc(dyn_m10B*sizeof(int *));  if (dyn_l10B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m10B; ++j) { dyn_l10B[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l10B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_10B = malloc(dyn_m10B*sizeof(int *));   if (dyn_sub_10B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10B[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m10B; ++j) { dyn_sub_10B[j] = malloc(3*sizeof(int));    if (dyn_sub_10B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_10B = malloc(m10B*sizeof(int));  if (dyn_up_10B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_10B[] malloc out of memory\n"); Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m10B; ++i) {
            for (j=0;j<10;++j) dyn_hc10B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l10B[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<3;++j) dyn_sub_10B[i][j]=-1;
        }
        for (i=0; i<m10B; ++i) {
            if (doSubClusts==1) dyn_up_10B[i]=-1;
        }
    }
    if (dyn_m10K!=-1) {
        dyn_hc10K = malloc(dyn_m10K*sizeof(int *)); if (dyn_hc10K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m10K; ++j) { dyn_hc10K[j] = malloc(10*sizeof(int)); if (dyn_hc10K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l10K = malloc(dyn_m10K*sizeof(int *));  if (dyn_l10K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m10K; ++j) { dyn_l10K[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l10K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m9K!=-1) {
            dyn_sub_10K = malloc(dyn_m10K*sizeof(int *));   if (dyn_sub_10K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10K[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m10K; ++j) { dyn_sub_10K[j] = malloc(1*sizeof(int));    if (dyn_sub_10K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m10K; ++i) {
            for (j=0;j<10;++j) dyn_hc10K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l10K[i][j]=-1;
            if (doSubClusts==1 && dyn_m9K!=-1) for (j=0;j<1;++j) dyn_sub_10K[i][j]=-1;
        }
    }
    if (dyn_m10W!=-1) {
        dyn_hc10W = malloc(dyn_m10W*sizeof(int *)); if (dyn_hc10W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10W[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m10W; ++j) { dyn_hc10W[j] = malloc(10*sizeof(int)); if (dyn_hc10W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc10W[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l10W = malloc(dyn_m10W*sizeof(int *));  if (dyn_l10W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10W[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m10W; ++j) { dyn_l10W[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l10W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l10W[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5b!=-1) {
            dyn_sub_10W = malloc(dyn_m10W*sizeof(int *));   if (dyn_sub_10W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10W[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m10W; ++j) { dyn_sub_10W[j] = malloc(6*sizeof(int));    if (dyn_sub_10W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_10W[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m10W; ++i) {
            for (j=0;j<10;++j) dyn_hc10W[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l10W[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5b!=-1) for (j=0;j<6;++j) dyn_sub_10W[i][j]=-1;
        }
    }
    if (dyn_m11A!=-1) {
        dyn_hc11A = malloc(dyn_m11A*sizeof(int *)); if (dyn_hc11A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11A; ++j) { dyn_hc11A[j] = malloc(11*sizeof(int)); if (dyn_hc11A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11A = malloc(dyn_m11A*sizeof(int *));  if (dyn_l11A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11A; ++j) { dyn_l11A[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m6A!=-1) {
            dyn_sub_11A = malloc(dyn_m11A*sizeof(int *));   if (dyn_sub_11A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11A[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11A; ++j) { dyn_sub_11A[j] = malloc(2*sizeof(int));    if (dyn_sub_11A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_11A = malloc(m11A*sizeof(int));  if (dyn_up_11A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_11A[] malloc out of memory\n"); Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m11A; ++i) {
            for (j=0;j<11;++j) dyn_hc11A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11A[i][j]=-1;
            if (doSubClusts==1 && dyn_m6A!=-1) for (j=0;j<2;++j) dyn_sub_11A[i][j]=-1;
        }
        for (i=0; i<m11A; ++i) {
            if (doSubClusts==1) dyn_up_11A[i]=-1;
        }
    }
    if (dyn_m11B!=-1) {
        dyn_hc11B = malloc(dyn_m11B*sizeof(int *)); if (dyn_hc11B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11B; ++j) { dyn_hc11B[j] = malloc(11*sizeof(int)); if (dyn_hc11B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11B = malloc(dyn_m11B*sizeof(int *));  if (dyn_l11B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11B; ++j) { dyn_l11B[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m9B!=-1) {
            dyn_sub_11B = malloc(dyn_m11B*sizeof(int *));   if (dyn_sub_11B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11B[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11B; ++j) { dyn_sub_11B[j] = malloc(1*sizeof(int));    if (dyn_sub_11B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m11B; ++i) {
            for (j=0;j<11;++j) dyn_hc11B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11B[i][j]=-1;
            if (doSubClusts==1 && dyn_m9B!=-1) for (j=0;j<1;++j) dyn_sub_11B[i][j]=-1;
        }
    }
    if (dyn_m11C!=-1) {
        dyn_hc11C = malloc(dyn_m11C*sizeof(int *)); if (dyn_hc11C==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11C[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11C; ++j) { dyn_hc11C[j] = malloc(11*sizeof(int)); if (dyn_hc11C[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11C[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11C = malloc(dyn_m11C*sizeof(int *));  if (dyn_l11C==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11C[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11C; ++j) { dyn_l11C[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11C[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11C[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_11C = malloc(dyn_m11C*sizeof(int *));   if (dyn_sub_11C==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11C[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11C; ++j) { dyn_sub_11C[j] = malloc(2*sizeof(int));    if (dyn_sub_11C[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11C[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_11C = malloc(m11C*sizeof(int));  if (dyn_up_11C==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_11C[] malloc out of memory\n"); Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m11C; ++i) {
            for (j=0;j<11;++j) dyn_hc11C[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11C[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<2;++j) dyn_sub_11C[i][j]=-1;
        }
        for (i=0; i<m11C; ++i) {
            if (doSubClusts==1) dyn_up_11C[i]=-1;
        }
    }
    if (dyn_m11E!=-1) {
        dyn_hc11E = malloc(dyn_m11E*sizeof(int *)); if (dyn_hc11E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11E[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11E; ++j) { dyn_hc11E[j] = malloc(11*sizeof(int)); if (dyn_hc11E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11E[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11E = malloc(dyn_m11E*sizeof(int *));  if (dyn_l11E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11E[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11E; ++j) { dyn_l11E[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11E[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_11E = malloc(dyn_m11E*sizeof(int *));   if (dyn_sub_11E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11E[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11E; ++j) { dyn_sub_11E[j] = malloc(3*sizeof(int));    if (dyn_sub_11E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11E[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m11E; ++i) {
            for (j=0;j<11;++j) dyn_hc11E[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11E[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<3;++j) dyn_sub_11E[i][j]=-1;
        }
    }
    if (dyn_m11F!=-1) {
        dyn_hc11F = malloc(dyn_m11F*sizeof(int *)); if (dyn_hc11F==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11F[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11F; ++j) { dyn_hc11F[j] = malloc(11*sizeof(int)); if (dyn_hc11F[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11F[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11F = malloc(dyn_m11F*sizeof(int *));  if (dyn_l11F==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11F[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11F; ++j) { dyn_l11F[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11F[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11F[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m6A!=-1) {
            dyn_sub_11F = malloc(dyn_m11F*sizeof(int *));   if (dyn_sub_11F==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11F[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11F; ++j) { dyn_sub_11F[j] = malloc(4*sizeof(int));    if (dyn_sub_11F[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11F[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        if (doSubClusts==1) {
            dyn_up_11F = malloc(m11F*sizeof(int));  if (dyn_up_11F==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_up_11F[] malloc out of memory\n"); Error_no_free(errMsg); }
        }
        for (i=0; i<dyn_m11F; ++i) {
            for (j=0;j<11;++j) dyn_hc11F[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11F[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m6A!=-1) for (j=0;j<4;++j) dyn_sub_11F[i][j]=-1;
        }
        for (i=0; i<m11F; ++i) {
            if (doSubClusts==1) dyn_up_11F[i]=-1;
        }
    }
    if (dyn_m11W!=-1) {
        dyn_hc11W = malloc(dyn_m11W*sizeof(int *)); if (dyn_hc11W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11W[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m11W; ++j) { dyn_hc11W[j] = malloc(11*sizeof(int)); if (dyn_hc11W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc11W[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l11W = malloc(dyn_m11W*sizeof(int *));  if (dyn_l11W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11W[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m11W; ++j) { dyn_l11W[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l11W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l11W[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m10B!=-1) {
            dyn_sub_11W = malloc(dyn_m11W*sizeof(int *));   if (dyn_sub_11W==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11W[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m11W; ++j) { dyn_sub_11W[j] = malloc(1*sizeof(int));    if (dyn_sub_11W[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_11W[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m11W; ++i) {
            for (j=0;j<11;++j) dyn_hc11W[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l11W[i][j]=-1;
            if (doSubClusts==1 && dyn_m10B!=-1) for (j=0;j<1;++j) dyn_sub_11W[i][j]=-1;
        }
    }
    if (dyn_m12A!=-1) {
        dyn_hc12A = malloc(dyn_m12A*sizeof(int *)); if (dyn_hc12A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m12A; ++j) { dyn_hc12A[j] = malloc(12*sizeof(int)); if (dyn_hc12A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l12A = malloc(dyn_m12A*sizeof(int *));  if (dyn_l12A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m12A; ++j) { dyn_l12A[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l12A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m11C!=-1) {
            dyn_sub_12A = malloc(dyn_m12A*sizeof(int *));   if (dyn_sub_12A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12A[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m12A; ++j) { dyn_sub_12A[j] = malloc(1*sizeof(int));    if (dyn_sub_12A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m12A; ++i) {
            for (j=0;j<12;++j) dyn_hc12A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l12A[i][j]=-1;
            if (doSubClusts==1 && dyn_m11C!=-1) for (j=0;j<1;++j) dyn_sub_12A[i][j]=-1;
        }
    }
    if (dyn_m12B!=-1) {
        dyn_hc12B = malloc(dyn_m12B*sizeof(int *)); if (dyn_hc12B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m12B; ++j) { dyn_hc12B[j] = malloc(12*sizeof(int)); if (dyn_hc12B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l12B = malloc(dyn_m12B*sizeof(int *));  if (dyn_l12B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m12B; ++j) { dyn_l12B[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l12B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_12B = malloc(dyn_m12B*sizeof(int *));   if (dyn_sub_12B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12B[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m12B; ++j) { dyn_sub_12B[j] = malloc(6*sizeof(int));    if (dyn_sub_12B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m12B; ++i) {
            for (j=0;j<12;++j) dyn_hc12B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l12B[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<6;++j) dyn_sub_12B[i][j]=-1;
        }
    }
    if (dyn_m12D!=-1) {
        dyn_hc12D = malloc(dyn_m12D*sizeof(int *)); if (dyn_hc12D==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12D[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m12D; ++j) { dyn_hc12D[j] = malloc(12*sizeof(int)); if (dyn_hc12D[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12D[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l12D = malloc(dyn_m12D*sizeof(int *));  if (dyn_l12D==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12D[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m12D; ++j) { dyn_l12D[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l12D[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12D[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_12D = malloc(dyn_m12D*sizeof(int *));   if (dyn_sub_12D==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12D[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m12D; ++j) { dyn_sub_12D[j] = malloc(4*sizeof(int));    if (dyn_sub_12D[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12D[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m12D; ++i) {
            for (j=0;j<12;++j) dyn_hc12D[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l12D[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<4;++j) dyn_sub_12D[i][j]=-1;
        }
    }
    if (dyn_m12E!=-1) {
        dyn_hc12E = malloc(dyn_m12E*sizeof(int *)); if (dyn_hc12E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12E[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m12E; ++j) { dyn_hc12E[j] = malloc(12*sizeof(int)); if (dyn_hc12E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12E[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l12E = malloc(dyn_m12E*sizeof(int *));  if (dyn_l12E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12E[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m12E; ++j) { dyn_l12E[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l12E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12E[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            dyn_sub_12E = malloc(dyn_m12E*sizeof(int *));   if (dyn_sub_12E==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12E[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m12E; ++j) { dyn_sub_12E[j] = malloc(3*sizeof(int));    if (dyn_sub_12E[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12E[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m12E; ++i) {
            for (j=0;j<12;++j) dyn_hc12E[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l12E[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1) for (j=0;j<3;++j) dyn_sub_12E[i][j]=-1;
        }
    }
    if (dyn_m12K!=-1) {
        dyn_hc12K = malloc(dyn_m12K*sizeof(int *)); if (dyn_hc12K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m12K; ++j) { dyn_hc12K[j] = malloc(12*sizeof(int)); if (dyn_hc12K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc12K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l12K = malloc(dyn_m12K*sizeof(int *));  if (dyn_l12K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m12K; ++j) { dyn_l12K[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l12K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l12K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m11A!=-1) {
            dyn_sub_12K = malloc(dyn_m12K*sizeof(int *));   if (dyn_sub_12K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12K[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m12K; ++j) { dyn_sub_12K[j] = malloc(1*sizeof(int));    if (dyn_sub_12K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_12K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m12K; ++i) {
            for (j=0;j<12;++j) dyn_hc12K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l12K[i][j]=-1;
            if (doSubClusts==1 && dyn_m11A!=-1) for (j=0;j<1;++j) dyn_sub_12K[i][j]=-1;
        }
    }
    if (dyn_m13A!=-1) {
        dyn_hc13A = malloc(dyn_m13A*sizeof(int *)); if (dyn_hc13A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13A[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m13A; ++j) { dyn_hc13A[j] = malloc(13*sizeof(int)); if (dyn_hc13A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13A[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l13A = malloc(dyn_m13A*sizeof(int *));  if (dyn_l13A==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13A[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m13A; ++j) { dyn_l13A[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l13A[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13A[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        for (i=0; i<dyn_m13A; ++i) {
            for (j=0;j<13;++j) dyn_hc13A[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l13A[i][j]=-1;
        }
    }
    if (dyn_m13B!=-1) {
        dyn_hc13B = malloc(dyn_m13B*sizeof(int *)); if (dyn_hc13B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13B[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m13B; ++j) { dyn_hc13B[j] = malloc(13*sizeof(int)); if (dyn_hc13B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l13B = malloc(dyn_m13B*sizeof(int *));  if (dyn_l13B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13B[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m13B; ++j) { dyn_l13B[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l13B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13B[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            dyn_sub_13B = malloc(dyn_m13B*sizeof(int *));   if (dyn_sub_13B==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_13B[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m13B; ++j) { dyn_sub_13B[j] = malloc(2*sizeof(int));    if (dyn_sub_13B[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_13B[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m13B; ++i) {
            for (j=0;j<13;++j) dyn_hc13B[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l13B[i][j]=-1;
            if (doSubClusts==1 && dyn_msp5c!=-1) for (j=0;j<2;++j) dyn_sub_13B[i][j]=-1;
        }
    }
    if (dyn_m13K!=-1) {
        dyn_hc13K = malloc(dyn_m13K*sizeof(int *)); if (dyn_hc13K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13K[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_m13K; ++j) { dyn_hc13K[j] = malloc(13*sizeof(int)); if (dyn_hc13K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hc13K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_l13K = malloc(dyn_m13K*sizeof(int *));  if (dyn_l13K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13K[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_m13K; ++j) { dyn_l13K[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_l13K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_l13K[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m11F!=-1) {
            dyn_sub_13K = malloc(dyn_m13K*sizeof(int *));   if (dyn_sub_13K==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_13K[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_m13K; ++j) { dyn_sub_13K[j] = malloc(3*sizeof(int));    if (dyn_sub_13K[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_13K[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_m13K; ++i) {
            for (j=0;j<13;++j) dyn_hc13K[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_l13K[i][j]=-1;
            if (doSubClusts==1 && dyn_m11F!=-1 && dyn_msp3c!=-1) for (j=0;j<3;++j) dyn_sub_13K[i][j]=-1;
        }
    }
    if (dyn_mFCC!=-1) {
        dyn_hcFCC = malloc(dyn_mFCC*sizeof(int *)); if (dyn_hcFCC==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcFCC[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_mFCC; ++j) { dyn_hcFCC[j] = malloc(13*sizeof(int)); if (dyn_hcFCC[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcFCC[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lFCC = malloc(dyn_mFCC*sizeof(int *));  if (dyn_lFCC==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lFCC[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_mFCC; ++j) { dyn_lFCC[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lFCC[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lFCC[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3b!=-1 && dyn_msp3c!=-1) {
            dyn_sub_FCC = malloc(dyn_mFCC*sizeof(int *));   if (dyn_sub_FCC==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_FCC[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_mFCC; ++j) { dyn_sub_FCC[j] = malloc(5*sizeof(int));    if (dyn_sub_FCC[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_FCC[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_mFCC; ++i) {
            for (j=0;j<13;++j) dyn_hcFCC[i][j]=-1; 
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lFCC[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3b!=-1 && dyn_msp3c!=-1) for (j=0;j<5;++j) dyn_sub_FCC[i][j]=-1;
        }
    }
    if (dyn_mHCP!=-1) {
        dyn_hcHCP = malloc(dyn_mHCP*sizeof(int *)); if (dyn_hcHCP==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcHCP[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_mHCP; ++j) { dyn_hcHCP[j] = malloc(13*sizeof(int)); if (dyn_hcHCP[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcHCP[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lHCP = malloc(dyn_mHCP*sizeof(int *));  if (dyn_lHCP==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lHCP[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_mHCP; ++j) { dyn_lHCP[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lHCP[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lHCP[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            dyn_sub_HCP = malloc(dyn_mHCP*sizeof(int *));   if (dyn_sub_HCP==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_HCP[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_mHCP; ++j) { dyn_sub_HCP[j] = malloc(3*sizeof(int));    if (dyn_sub_HCP[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_HCP[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_mHCP; ++i) {
            for (j=0;j<13;++j) dyn_hcHCP[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lHCP[i][j]=-1;
            if (doSubClusts==1 && dyn_msp3c!=-1) for (j=0;j<3;++j) dyn_sub_HCP[i][j]=-1;
        }
    }
    if (dyn_mBCC_9!=-1) {
        dyn_hcBCC_9 = malloc(dyn_mBCC_9*sizeof(int *)); if (dyn_hcBCC_9==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcBCC_9[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_mBCC_9; ++j) { dyn_hcBCC_9[j] = malloc(15*sizeof(int)); if (dyn_hcBCC_9[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcBCC_9[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        dyn_lBCC_9 = malloc(dyn_mBCC_9*sizeof(int *));  if (dyn_lBCC_9==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lBCC_9[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_mBCC_9; ++j) { dyn_lBCC_9[j] = malloc((2*initNoLifetimes+1)*sizeof(int));   if (dyn_lBCC_9[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lBCC_9[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_msp4b!=-1 && dyn_m6A!=-1) {
            dyn_sub_BCC_9 = malloc(dyn_mBCC_9*sizeof(int *));   if (dyn_sub_BCC_9==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_BCC_9[] malloc out of memory\n");   Error_no_free(errMsg); }
            for (j=0; j<dyn_mBCC_9; ++j) { dyn_sub_BCC_9[j] = malloc(6*sizeof(int));    if (dyn_sub_BCC_9[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_BCC_9[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_mBCC_9; ++i) {
            for (j=0; j<9; ++j) dyn_hcBCC_9[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lBCC_9[i][j]=-1;
            if (doSubClusts==1 && dyn_msp4b!=-1 && dyn_m6A!=-1) for (j=0;j<6;++j) dyn_sub_BCC_9[i][j]=-1;
        }
    }
    if (dyn_mBCC_15!=-1) {
        dyn_hcBCC_15 = malloc(dyn_mBCC_15*sizeof(int *));   if (dyn_hcBCC_15==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcBCC_15[] malloc out of memory\n"); Error_no_free(errMsg); }
        for (j=0; j<dyn_mBCC_15; ++j) { dyn_hcBCC_15[j] = malloc(15*sizeof(int));   if (dyn_hcBCC_15[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_hcBCC_15[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        dyn_lBCC_15 = malloc(dyn_mBCC_15*sizeof(int *));    if (dyn_lBCC_15==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lBCC_15[] malloc out of memory\n");   Error_no_free(errMsg); }
        for (j=0; j<dyn_mBCC_15; ++j) { dyn_lBCC_15[j] = malloc((2*initNoLifetimes+1)*sizeof(int)); if (dyn_lBCC_15[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_lBCC_15[][] malloc out of memory\n");  Error_no_free(errMsg); } }
        if (doSubClusts==1 && dyn_m6A!=-1) {
            dyn_sub_BCC_15 = malloc(dyn_mBCC_15*sizeof(int *)); if (dyn_sub_BCC_15==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_BCC_15[] malloc out of memory\n"); Error_no_free(errMsg); }
            for (j=0; j<dyn_mBCC_15; ++j) { dyn_sub_BCC_15[j] = malloc(6*sizeof(int));  if (dyn_sub_BCC_15[j]==NULL) { sprintf(errMsg,"Setup_InitDynamicVars(): dyn_sub_BCC_15[][] malloc out of memory\n");    Error_no_free(errMsg); } }
        }
        for (i=0; i<dyn_mBCC_15; ++i) {
            for (j=0; j<15; ++j) dyn_hcBCC_15[i][j]=-1;
            for (j=0;j<(2*initNoLifetimes+1);++j) dyn_lBCC_15[i][j]=-1;
            if (doSubClusts==1 && dyn_m6A!=-1) for (j=0;j<6;++j) dyn_sub_BCC_15[i][j]=-1;
        }
    }
}

void Setup_FreeDynamicVars() {  // Free bond detection variables
    int j;
    int printdebugs=0;
    
    if (printdebugs==1) printf("hi1\n"); if (dyn_msp3!=-1) {
        for (j=0; j<dyn_msp3; ++j) { 
            free(dyn_sp3[j]);
            free(dyn_lsp3[j]);
        }
        free(dyn_sp3);
        free(dyn_lsp3);
    }
    if (printdebugs==1) printf("hi2\n"); if (dyn_msp3a!=-1) {
        for (j=0; j<dyn_msp3a; ++j) { 
            free(dyn_sp3a[j]);
            free(dyn_lsp3a[j]);
        }
        free(dyn_sp3a);
        free(dyn_lsp3a);
    }
    if (printdebugs==1) printf("hi3\n"); if (dyn_msp3b!=-1) {
        for (j=0; j<dyn_msp3b; ++j) { 
            free(dyn_sp3b[j]);
            free(dyn_lsp3b[j]);
        }
        free(dyn_sp3b);
        free(dyn_lsp3b);
        if (doSubClusts==1) free(dyn_up_sp3b);
    }
    if (printdebugs==1) printf("hi4\n"); if (dyn_msp3c!=-1) {
        for (j=0; j<dyn_msp3c; ++j) { 
            free(dyn_sp3c[j]);
            free(dyn_lsp3c[j]);
        }
        free(dyn_sp3c);
        free(dyn_lsp3c);
        if (doSubClusts==1) free(dyn_up_sp3c);
    }
    if (printdebugs==1) printf("hi5\n"); if (dyn_msp4!=-1) {
        for (j=0; j<dyn_msp4; ++j) { 
            free(dyn_sp4[j]);
            free(dyn_lsp4[j]);
        }
        free(dyn_sp4);
        free(dyn_lsp4);
    }
    if (printdebugs==1) printf("hi6\n"); if (dyn_msp4a!=-1) {
        for (j=0; j<dyn_msp4a; ++j) { 
            free(dyn_sp4a[j]);
            free(dyn_lsp4a[j]);
        }
        free(dyn_sp4a);
        free(dyn_lsp4a);
    }
    if (printdebugs==1) printf("hi7\n"); if (dyn_msp4b!=-1) {
        for (j=0; j<dyn_msp4b; ++j) { 
            free(dyn_sp4b[j]);
            free(dyn_lsp4b[j]);
        }
        free(dyn_sp4b);
        free(dyn_lsp4b);
        if (doSubClusts==1) free(dyn_up_sp4b);
    }
    if (printdebugs==1) printf("hi8\n"); if (dyn_m6A!=-1) {
        for (j=0; j<dyn_m6A; ++j) { 
            free(dyn_hc6A[j]);
            free(dyn_l6A[j]);
        }
        free(dyn_hc6A);
        free(dyn_l6A);
        if (doSubClusts==1) free(dyn_up_sp4c);
    }
    if (printdebugs==1) printf("hi9\n"); if (dyn_msp5!=-1) {
        if (printdebugs==1) printf("hi9.1\n"); for (j=0; j<dyn_msp5; ++j) { 
            if (printdebugs==1) printf("hi9.1.1\n"); free(dyn_sp5[j]);
            if (printdebugs==1) printf("hi9.1.2 %d\n",j); free(dyn_lsp5[j]); if (printdebugs==1) printf("hi9.1.3\n");
        }
        if (printdebugs==1) printf("hi9.2\n"); free(dyn_sp5);
        if (printdebugs==1) printf("hi9.3\n"); free(dyn_lsp5); if (printdebugs==1) printf("hi9.4\n");
    }
    if (printdebugs==1) printf("hi10\n"); if (dyn_msp5a!=-1) {
        for (j=0; j<dyn_msp5a; ++j) { 
            free(dyn_sp5a[j]);
            free(dyn_lsp5a[j]);
        }
        free(dyn_sp5a);
        free(dyn_lsp5a);
    }
    if (printdebugs==1) printf("hi11\n"); if (dyn_msp5b!=-1) {
        for (j=0; j<dyn_msp5b; ++j) { 
            free(dyn_sp5b[j]);
            free(dyn_lsp5b[j]);
        }
        free(dyn_sp5b);
        free(dyn_lsp5b);
        if (doSubClusts==1) free(dyn_up_sp5b);
    }
    if (printdebugs==1) printf("hi12\n"); if (dyn_msp5c!=-1) {
        for (j=0; j<dyn_msp5c; ++j) { 
            free(dyn_sp5c[j]);
            free(dyn_lsp5c[j]);
        }
        free(dyn_sp5c);
        free(dyn_lsp5c);
        if (doSubClusts==1) free(dyn_up_sp5c);
    }
    if (printdebugs==1) printf("hi13\n"); if (dyn_m6Z!=-1) {
        for (j=0; j<dyn_m6Z; ++j) { 
            free(dyn_hc6Z[j]);
            free(dyn_l6Z[j]);
        }
        free(dyn_hc6Z);
        free(dyn_l6Z);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m6Z; ++j) free(dyn_sub_6Z[j]); 
            free(dyn_sub_6Z);
        }
    }
    if (printdebugs==1) printf("hi13.1\n"); if (dyn_m7K!=-1) {
        for (j=0; j<dyn_m7K; ++j) { 
            free(dyn_hc7K[j]);
            free(dyn_l7K[j]);
        }
        free(dyn_hc7K);
        free(dyn_l7K);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m7K; ++j) free(dyn_sub_7K[j]); 
            free(dyn_sub_7K);
        }
    }
    if (printdebugs==1) printf("hi14\n"); if (dyn_m8A!=-1) {
        for (j=0; j<dyn_m8A; ++j) { 
            free(dyn_hc8A[j]);
            free(dyn_l8A[j]);
        }
        free(dyn_hc8A);
        free(dyn_l8A);
        if (doSubClusts==1 && dyn_msp5b!=-1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m8A; ++j) free(dyn_sub_8A[j]); 
            free(dyn_sub_8A);
        }
    }
    if (printdebugs==1) printf("hi15\n"); if (dyn_m8B!=-1) {
        for (j=0; j<dyn_m8B; ++j) { 
            free(dyn_hc8B[j]);
            free(dyn_l8B[j]);
        }
        free(dyn_hc8B);
        free(dyn_l8B);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m8B; ++j) free(dyn_sub_8B[j]); 
            free(dyn_sub_8B);
        }
    }
    if (printdebugs==1) printf("hi15.5\n"); if (dyn_m8K!=-1) {
        for (j=0; j<dyn_m8K; ++j) { 
            free(dyn_hc8K[j]);
            free(dyn_l8K[j]);
        }
        free(dyn_hc8K);
        free(dyn_l8K);
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            for (j=0; j<dyn_m8K; ++j) free(dyn_sub_8K[j]); 
            free(dyn_sub_8K);
        }
    }
    if (printdebugs==1) printf("hi16\n"); if (dyn_m9A!=-1) {
        for (j=0; j<dyn_m9A; ++j) { 
            free(dyn_hc9A[j]);
            free(dyn_l9A[j]);
        }
        free(dyn_hc9A);
        free(dyn_l9A);
        if (doSubClusts==1 && dyn_msp4b!=-1) {
            for (j=0; j<dyn_m9A; ++j) free(dyn_sub_9A[j]); 
            free(dyn_sub_9A);
        }
    }
    if (printdebugs==1) printf("hi17\n"); if (dyn_m9B!=-1) {
        for (j=0; j<dyn_m9B; ++j) { 
            free(dyn_hc9B[j]);
            free(dyn_l9B[j]);
        }
        free(dyn_hc9B);
        free(dyn_l9B);
        if (doSubClusts==1) free(dyn_up_9B);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m9B; ++j) free(dyn_sub_9B[j]); 
            free(dyn_sub_9B);
        }
    }
    if (printdebugs==1) printf("hi17.5\n"); if (dyn_m9K!=-1) {
        for (j=0; j<dyn_m9K; ++j) { 
            free(dyn_hc9K[j]);
            free(dyn_l9K[j]);
        }
        free(dyn_hc9K);
        free(dyn_l9K);
        if (doSubClusts==1) free(dyn_up_9K);
        if (doSubClusts==1 && dyn_m6A!=-1) {
            for (j=0; j<dyn_m9K; ++j) free(dyn_sub_9K[j]); 
            free(dyn_sub_9K);
        }
    }
    if (printdebugs==1) printf("hi18\n"); if (dyn_m10A!=-1) {
        for (j=0; j<dyn_m10A; ++j) { 
            free(dyn_hc10A[j]);
            free(dyn_l10A[j]);
        }
        free(dyn_hc10A);
        free(dyn_l10A);
        if (doSubClusts==1 && dyn_msp4b!=-1) {
            for (j=0; j<dyn_m10A; ++j) free(dyn_sub_10A[j]);
            free(dyn_sub_10A);
        }
    }
    if (printdebugs==1) printf("hi19\n"); if (dyn_m10B!=-1) {
        for (j=0; j<dyn_m10B; ++j) { 
            free(dyn_hc10B[j]);
            free(dyn_l10B[j]);
        }
        free(dyn_hc10B);
        free(dyn_l10B);
        if (doSubClusts==1) free(dyn_up_10B);
        if (doSubClusts==1 && dyn_msp5c!=-1 && dyn_m9B!=-1) {
            for (j=0; j<dyn_m10B; ++j) free(dyn_sub_10B[j]); 
            free(dyn_sub_10B);
        }
    }
    if (printdebugs==1) printf("hi19.25\n"); if (dyn_m10K!=-1) {
        for (j=0; j<dyn_m10K; ++j) { 
            free(dyn_hc10K[j]);
            free(dyn_l10K[j]);
        }
        free(dyn_hc10K);
        free(dyn_l10K);
        if (doSubClusts==1 && dyn_m9K!=-1) {
            for (j=0; j<dyn_m10K; ++j) free(dyn_sub_10K[j]); 
            free(dyn_sub_10K);
        }
    }
    if (printdebugs==1) printf("hi19.5\n"); if (dyn_m10W!=-1) {
        for (j=0; j<dyn_m10W; ++j) { 
            free(dyn_hc10W[j]);
            free(dyn_l10W[j]);
        }
        free(dyn_hc10W);
        free(dyn_l10W);
        if (doSubClusts==1 && dyn_msp5b!=-1) {
            for (j=0; j<dyn_m10W; ++j) free(dyn_sub_10W[j]); 
            free(dyn_sub_10W);
        }
    }
    if (printdebugs==1) printf("hi20\n"); if (dyn_m11A!=-1) {
        for (j=0; j<dyn_m11A; ++j) { 
            free(dyn_hc11A[j]);
            free(dyn_l11A[j]);
        }
        free(dyn_hc11A);
        free(dyn_l11A);
        if (doSubClusts==1) free(dyn_up_11A);
        if (doSubClusts==1 && dyn_m6A!=-1) {
            for (j=0; j<dyn_m11A; ++j) free(dyn_sub_11A[j]); 
            free(dyn_sub_11A);
        }
    }
    if (printdebugs==1) printf("hi21\n"); if (dyn_m11B!=-1) {
        for (j=0; j<dyn_m11B; ++j) { 
            free(dyn_hc11B[j]);
            free(dyn_l11B[j]);
        }
        free(dyn_hc11B);
        free(dyn_l11B);
        if (doSubClusts==1 && dyn_m9B!=-1) {
            for (j=0; j<dyn_m11B; ++j) free(dyn_sub_11B[j]);
            free(dyn_sub_11B);
        }
    }
    if (printdebugs==1) printf("hi22\n"); if (dyn_m11C!=-1) {
        for (j=0; j<dyn_m11C; ++j) { 
            free(dyn_hc11C[j]);
            free(dyn_l11C[j]);
        }
        free(dyn_hc11C);
        free(dyn_l11C);
        if (doSubClusts==1) free(dyn_up_11C);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m11C; ++j) free(dyn_sub_11C[j]); 
            free(dyn_sub_11C);
        }
    }
    if (printdebugs==1) printf("hi23\n"); if (dyn_m11E!=-1) {
        for (j=0; j<dyn_m11E; ++j) { 
            free(dyn_hc11E[j]);
            free(dyn_l11E[j]);
        }
        free(dyn_hc11E);
        free(dyn_l11E);
        if (doSubClusts==1 && dyn_msp5c!=-1 && dyn_m9B!=-1) {
            for (j=0; j<dyn_m11E; ++j) free(dyn_sub_11E[j]); 
            free(dyn_sub_11E);
        }
    }
    if (printdebugs==1) printf("hi24\n"); if (dyn_m11F!=-1) {
        for (j=0; j<dyn_m11F; ++j) { 
            free(dyn_hc11F[j]);
            free(dyn_l11F[j]);
        }
        free(dyn_hc11F);
        free(dyn_l11F);
        if (doSubClusts==1) free(dyn_up_11F);
        if (doSubClusts==1 && dyn_msp3c!=-1 && dyn_m6A!=-1) {
            for (j=0; j<dyn_m11F; ++j) free(dyn_sub_11F[j]); 
            free(dyn_sub_11F);
        }
    }
    if (printdebugs==1) printf("hi24.2\n"); if (dyn_m11W!=-1) {
        for (j=0; j<dyn_m11W; ++j) { 
            free(dyn_hc11W[j]);
            free(dyn_l11W[j]);
        }
        free(dyn_hc11W);
        free(dyn_l11W);
        if (doSubClusts==1 && dyn_m10B!=-1) {
            for (j=0; j<dyn_m11W; ++j) free(dyn_sub_11W[j]); 
            free(dyn_sub_11W);
        }
    }
    if (printdebugs==1) printf("hi25\n"); if (dyn_m12A!=-1) {
        for (j=0; j<dyn_m12A; ++j) { 
            free(dyn_hc12A[j]);
            free(dyn_l12A[j]);
        }
        free(dyn_hc12A);
        free(dyn_l12A);
        if (doSubClusts==1 && dyn_m11C!=-1) {
            for (j=0; j<dyn_m12A; ++j) free(dyn_sub_12A[j]);
            free(dyn_sub_12A);
        }
    }
    if (printdebugs==1) printf("hi26\n"); if (dyn_m12B!=-1) {
        for (j=0; j<dyn_m12B; ++j) { 
            free(dyn_hc12B[j]);
            free(dyn_l12B[j]);
        }
        free(dyn_hc12B);
        free(dyn_l12B);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m12B; ++j) free(dyn_sub_12B[j]); 
            free(dyn_sub_12B);
        }
    }
    if (printdebugs==1) printf("hi27\n"); if (dyn_m12D!=-1) {
        for (j=0; j<dyn_m12D; ++j) { 
            free(dyn_hc12D[j]);
            free(dyn_l12D[j]);
        }
        free(dyn_hc12D);
        free(dyn_l12D);
        if (doSubClusts==1 && dyn_msp5c!=-1 && dyn_m11E!=-1) {
            for (j=0; j<dyn_m12D; ++j) free(dyn_sub_12D[j]); 
            free(dyn_sub_12D);
        }
    }
    if (printdebugs==1) printf("hi28\n"); if (dyn_m12E!=-1) {
        for (j=0; j<dyn_m12E; ++j) { 
            free(dyn_hc12E[j]);
            free(dyn_l12E[j]);
        }
        free(dyn_hc12E);
        free(dyn_l12E);
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            for (j=0; j<dyn_m12E; ++j) free(dyn_sub_12E[j]);
            free(dyn_sub_12E);
        }
    }
    if (printdebugs==1) printf("hi28.5\n"); if (dyn_m12K!=-1) {
        for (j=0; j<dyn_m12K; ++j) { 
            free(dyn_hc12K[j]);
            free(dyn_l12K[j]);
        }
        free(dyn_hc12K);
        free(dyn_l12K);
        if (doSubClusts==1 && dyn_m11A!=-1) {
            for (j=0; j<dyn_m12K; ++j) free(dyn_sub_12K[j]);
            free(dyn_sub_12K);
        }
    }
    if (printdebugs==1) printf("hi29\n"); if (dyn_m13A!=-1) {
        for (j=0; j<dyn_m13A; ++j) {
            free(dyn_hc13A[j]);
            free(dyn_l13A[j]);
        }
        free(dyn_hc13A);
        free(dyn_l13A);
        
    }
    if (printdebugs==1) printf("hi30\n"); if (dyn_m13B!=-1) {
        for (j=0; j<dyn_m13B; ++j) { 
            free(dyn_hc13B[j]);
            free(dyn_l13B[j]);
        }
        free(dyn_hc13B);
        free(dyn_l13B);
        if (doSubClusts==1 && dyn_msp5c!=-1) {
            for (j=0; j<dyn_m13B; ++j) free(dyn_sub_13B[j]); 
            free(dyn_sub_13B);
        }
    }
    if (printdebugs==1) printf("hi28\n"); if (dyn_m13K!=-1) {
        for (j=0; j<dyn_m13K; ++j) { 
            free(dyn_hc13K[j]);
            free(dyn_l13K[j]);
        }
        free(dyn_hc13K);
        free(dyn_l13K);
        if (doSubClusts==1 && dyn_m11F!=-1) {
            for (j=0; j<dyn_m13K; ++j) free(dyn_sub_13K[j]);
            free(dyn_sub_13K);
        }
    }
    if (printdebugs==1) printf("hi31\n"); if (dyn_mFCC!=-1) {
        for (j=0; j<dyn_mFCC; ++j) { 
            free(dyn_hcFCC[j]);
            free(dyn_lFCC[j]);
        }
        free(dyn_hcFCC);
        free(dyn_lFCC);
        if (doSubClusts==1 && dyn_msp3b!=-1 && dyn_msp3c!=-1) {
            for (j=0; j<dyn_mFCC; ++j) free(dyn_sub_FCC[j]); 
            free(dyn_sub_FCC);
        }
    }
    if (printdebugs==1) printf("hi32\n"); if (dyn_mHCP!=-1) {
        for (j=0; j<dyn_mHCP; ++j) { 
            free(dyn_hcHCP[j]);
            free(dyn_lHCP[j]);
        }
        free(dyn_hcHCP);
        free(dyn_lHCP);
        if (doSubClusts==1 && dyn_msp3c!=-1) {
            for (j=0; j<dyn_mHCP; ++j) free(dyn_sub_HCP[j]); 
            free(dyn_sub_HCP);
        }
    }
    if (printdebugs==1) printf("hi33\n"); if (dyn_mBCC_9!=-1) {
        for (j=0; j<dyn_mBCC_9; ++j) { 
            free(dyn_hcBCC_9[j]);
            free(dyn_lBCC_9[j]);
        }
        free(dyn_hcBCC_9);
        free(dyn_lBCC_9);
        if (doSubClusts==1 && dyn_m6A!=-1) {
            for (j=0; j<dyn_mBCC_9; ++j) free(dyn_sub_BCC_9[j]); 
            free(dyn_sub_BCC_9);
        }
    }
    if (printdebugs==1) printf("hi34\n"); if (dyn_mBCC_15!=-1) {
        for (j=0; j<dyn_mBCC_15; ++j) { 
            free(dyn_hcBCC_15[j]);
            free(dyn_lBCC_15[j]);
        }
        free(dyn_hcBCC_15);
        free(dyn_lBCC_15);
        if (doSubClusts==1 && dyn_m6A!=-1) {
            for (j=0; j<dyn_mBCC_15; ++j) free(dyn_sub_BCC_15[j]); 
            free(dyn_sub_BCC_15);
        }
    }
    if (printdebugs==1) printf("hi35\n"); 
}