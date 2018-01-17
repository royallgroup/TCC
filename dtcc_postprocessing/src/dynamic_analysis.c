/* TCC analysis */
/* Alex Malins, Bristol Centre for Complexity Sciences, University of Bristol, United Kingdom */
/* 2012 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dynamic_analysis.h"
#include "iniparser.h"

void Error(char *msg) {
	printf("\n%s\n",msg);
	exit(1);
} // Quit

//// START: Setup routines
int Setup_GetFirstIntFromLine(FILE *stream) { // returns first integer from line in file stream
	char input[10000], errMsg[1000];
	char *pch;
	if (fgets(input,10000,stream)==NULL) {
		sprintf(errMsg,"Setup_GetFirstIntFromLine(): end of input file reached\n");
		Error(errMsg);
	}
	pch=strtok(input," ");

	return atoi(pch);
}

double Setup_GetFirstDoubleFromLine(FILE *stream) { // returns first double from line in file stream
	char input[10000], errMsg[1000];
	char *pch;
	if (fgets(input,10000,stream)==NULL) {
		sprintf(errMsg,"Setup_GetFirstDoubleFromLine(): end of input file reached\n");
		Error(errMsg);
	}
	pch=strtok(input," ");
	return atof(pch);
}

void Setup_ReadIniFile(char *filename) { // reads analysis code params file

    dictionary  *   ini ;
    char errMsg[1000];
    double RHO;
    int i;

    fXmolName=malloc(500*sizeof(char)); if (fXmolName==NULL) { sprintf(errMsg,"Setup_InitStaticVars(): fXmolName[] malloc out of memory\n");   Error(errMsg); }

	ini = iniparser_load(filename);
    if (ini==NULL) {
        sprintf(errMsg,"Setup_ReadIniFile(): Error opening file %s",filename);  // Always test file open
        Error(errMsg);
    }

	printf("\n%s reading inputfile:\n",filename);

    //box
    ISNOTCUBIC = iniparser_getint(ini, "box:box_type", -1);
    strcpy(fBoxSizeName, (char*)iniparser_getstring(ini, "box:box_name", "-1"));    
    

	//run
    strcpy(fXmolName, (char*)iniparser_getstring(ini, "run:xyzfilename", "-1"));
    FRAMES = iniparser_getint(ini, "run:frames", -1);
    N = iniparser_getint(ini, "run:num_particles", -1);
    NA = iniparser_getint(ini, "run:numA_particles", -1);
    STARTFROM = iniparser_getint(ini, "run:start_from", -1);
    SAMPLEFREQ = iniparser_getint(ini, "run:sample_freqency", -1);
    RHO = iniparser_getdouble(ini, "run:number_density", -1);
    TSTART = iniparser_getdouble(ini, "run:simulationstarttime", -1);
    FRAMETSTEP = iniparser_getdouble(ini, "run:simulationtimestep", -1);
    TFINAL = iniparser_getdouble(ini, "run:simulationendtime", -1);


    //simulation
    rcutAA = iniparser_getdouble(ini, "simulation:rcutAA", -1);
    rcutAB = iniparser_getdouble(ini, "simulation:rcutAB", -1);
    rcutBB = iniparser_getdouble(ini, "simulation:rcutBB", -1);
    Vor = iniparser_getboolean(ini, "simulation:bond_type", -1);
    PBCs = iniparser_getboolean(ini, "simulation:pbcs", -1);
    fc = iniparser_getdouble(ini, "simulation:voronoi_parameter", -1);
    nB = iniparser_getint(ini, "simulation:num_bonds", -1);

    //output
    binWidth = iniparser_getdouble(ini, "output:bin_width", -1);
    doSubClusts = iniparser_getboolean(ini, "output:subclusters", -1);

	//extra
    talpha = iniparser_getdouble(ini, "extra:alpha_time", -1);
    PRINTINFO = iniparser_getboolean(ini, "extra:debug", -1);

    // DTCC post processing
    doDynamicAnalysis = iniparser_getboolean(ini, "dynamicpp:dynamic_clusters", -1);
	go_away_for	= iniparser_getint(ini, "dynamicpp:max_dissapear", -1);
	doSplits = iniparser_getboolean(ini, "dynamicpp:remove_unbonded", -1);
	
	doDynFileRewrite = iniparser_getboolean(ini, "dynamicpp:rewrite_dyn", -1);
	doWriteNewClustsFile = iniparser_getboolean(ini, "dynamicpp:write_new", -1);
	steady_Nclus = iniparser_getint(ini, "dynamicpp:steady_state", -1);
	doWriteLifetimeDistro = iniparser_getboolean(ini, "dynamicpp:lifetimes", -1);
	doWriteIn = iniparser_getboolean(ini, "dynamicpp:raw_trajectory", -1);
	t_h = iniparser_getint(ini, "dynamicpp:frames_displacement", -1);
	doDisplacementDistro = iniparser_getboolean(ini, "dynamicpp:displacement_dist", -1);
	doWriteCtsIn = iniparser_getboolean(ini, "dynamicpp:raw_trajclusters", -1);
	doWriteMSDClusNonClusCts = iniparser_getboolean(ini, "dynamicpp:msd_clusters", -1);
	max_regions = iniparser_getint(ini, "dynamicpp:clusterregions", -1);
	doRegions = iniparser_getboolean(ini, "dynamicpp:region_size", -1);

    //Post reading calculations
	if (doDynamicAnalysis==1) {
		start_from=0+go_away_for+1;
		end_by=FRAMES-1-go_away_for-1;
		useable_frames=0;
		for (i=0; i<FRAMES; i++) {
			if (i<start_from || i>end_by) continue;
			useable_frames++;
		}
        if (useable_frames < 1) {
            sprintf(errMsg,"Setup_ReadInputParams(): There is less than 1 usable frame. Increase number of frames to analyse or decreases steady state frames.");
            Error(errMsg);
        }

		percent_start=0+steady_Nclus;
		percent_end=FRAMES-1-steady_Nclus;
		useable_percent=percent_end-percent_start+1;
	}
	else {
		start_from=percent_start=0;
		end_by=percent_end=FRAMES-1;
		useable_frames=FRAMES;
	}

	if (NA<N) doBinary=1;
	else doBinary=0;

	if (ISNOTCUBIC==0) {
		sidex=sidey=sidez=pow((double)N/RHO, 1.0/3.0);
		halfSidex=halfSidey=halfSidez=sidex/2.0;
		printf("Cubic box side length = %.15lg\n\n",sidex);
	}

	//Print out values
    printf("frames to sample from xmol %d\n",FRAMES);
    printf("start sampling from frame %d in xmol file\n",STARTFROM);
    printf("sampling frequency of frames %d\n",SAMPLEFREQ);
    printf("Usable frames = %d\n", useable_frames);
	printf("only consider clusters which are created on or after %d and end on or before %d of FRAMES %d\n",start_from,end_by,FRAMES);
	printf("therefore useable_frames is %d\n",useable_frames);
	printf("calculate Nclus(tau_l)/N between frames %d and %d\n",percent_start,percent_end);
	printf("therefore useable frames for Nclus(tau_l)/N is %d\n",useable_percent);
	if (doDynamicAnalysis==1) printf("doDynamicAnalysis %d therefore removing disappearences of %d or fewer frames\n",doDynamicAnalysis,go_away_for);
	else if (doDynamicAnalysis==0)  printf("doDynamicAnalysis %d therefore cluster disappearences are not removed\n",doDynamicAnalysis);
	if (doDynamicAnalysis==0) doSplits=0;
	printf("remove instances where subcluster becomes unbonded from cluster during a disappearence event? %d\n",doSplits);
	printf("rewrite dyn_... files for debugging? %d\n",doDynFileRewrite);
	printf("write processed clusters to new_... file %d\n",doWriteNewClustsFile);
	printf("number of frames for N_clus/N to reach steady state - needs to be of order of length of longest lived cluster (frames) %d\n",steady_Nclus);
	printf("output lifetime distributions to lives_... file %d\n",doWriteLifetimeDistro);
	printf("write raw style trajectory of particles in clusters %d\n",doWriteIn);
	printf("t_h - number of frames to write displacement distributions, particles continuously in clusters, etc %d\n",t_h);
	printf("write displacement distribution of clustered particles from t to t+%d frames %d\n",t_h,doDisplacementDistro);
	printf("write raw style trajectory of particles in continuously in clusters over t_h %d - %d\n",t_h,doWriteCtsIn);
	printf("write MSD of particles in continuously in clusters over t - %d\n",doWriteMSDClusNonClusCts);
	printf("maximum number of regions of any cluster type %d\n",max_regions);
	printf("average size of regions of clusters, i.e. look at domain properties of regions of clusters %d\n",doRegions);


	iniparser_freedict(ini);

}


void Setup_ReadBox(FILE *readIn)  {
    int sweep;
    double tilt;

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

void Setup_ReadDynamicDat(char *filename) { // Initialize lots of important variables
	char errMsg[1000];
	FILE *fin;

	printf("reading dynamics memory parameters from %s\n",filename);
	fin=fopen(filename,"r");
	if (fin==NULL)  {
		sprintf(errMsg,"Setup_ReadDynamicDat() : Error opening file %s",filename);	// Always test file open
		Error(errMsg);
	}

	dyn_msp3=Setup_GetFirstIntFromLine(fin);
	dyn_msp3a=Setup_GetFirstIntFromLine(fin);
	dyn_msp3b=Setup_GetFirstIntFromLine(fin);
	dyn_msp3c=Setup_GetFirstIntFromLine(fin);
	printf("dyn_msp3 %d dyn_msp3a %d dyn_msp3b %d dyn_msp3c %d\n",dyn_msp3,dyn_msp3a,dyn_msp3b,dyn_msp3c);
	dyn_msp4=Setup_GetFirstIntFromLine(fin);
	dyn_msp4a=Setup_GetFirstIntFromLine(fin);
	dyn_msp4b=Setup_GetFirstIntFromLine(fin);
	dyn_m6A=Setup_GetFirstIntFromLine(fin);
	printf("dyn_msp4 %d dyn_msp4a %d dyn_msp4b %d dyn_m6A %d\n",dyn_msp4,dyn_msp4a,dyn_msp4b,dyn_m6A);
	dyn_msp5=Setup_GetFirstIntFromLine(fin);
	dyn_msp5a=Setup_GetFirstIntFromLine(fin);
	dyn_msp5b=Setup_GetFirstIntFromLine(fin);
	dyn_msp5c=Setup_GetFirstIntFromLine(fin);
	printf("dyn_msp5 %d dyn_msp5a %d dyn_msp5b %d dyn_msp5c %d\n",dyn_msp5,dyn_msp5a,dyn_msp5b,dyn_msp5c);
	dyn_m6Z=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m6Z %d\n",dyn_m6Z);
	dyn_m7K=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m7K %d\n",dyn_m7K);
	dyn_m8A=Setup_GetFirstIntFromLine(fin);
	dyn_m8B=Setup_GetFirstIntFromLine(fin);
	dyn_m8K=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m8A %d dyn_m8B %d dyn_m8K %d\n",dyn_m8A,dyn_m8B,dyn_m8K);
	dyn_m9A=Setup_GetFirstIntFromLine(fin);
	dyn_m9B=Setup_GetFirstIntFromLine(fin);
	dyn_m9K=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m9A %d dyn_m9B %d dyn_m9K %d\n",dyn_m9A,dyn_m9B,dyn_m9K);
	dyn_m10A=Setup_GetFirstIntFromLine(fin);
	dyn_m10B=Setup_GetFirstIntFromLine(fin);
	dyn_m10K=Setup_GetFirstIntFromLine(fin);
	dyn_m10W=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m10A %d dyn_m10B %d dyn_m10K %d dyn_m10W %d\n",dyn_m10A,dyn_m10B,dyn_m10K,dyn_m10W);
	dyn_m11A=Setup_GetFirstIntFromLine(fin);
	dyn_m11B=Setup_GetFirstIntFromLine(fin);
	dyn_m11C=Setup_GetFirstIntFromLine(fin);
	dyn_m11E=Setup_GetFirstIntFromLine(fin);
	dyn_m11F=Setup_GetFirstIntFromLine(fin);
	dyn_m11W=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m11A %d dyn_m11B %d dyn_m11C %d dyn_m11E %d dyn_m11F %d dyn_m11W %d\n",dyn_m11A,dyn_m11B,dyn_m11C,dyn_m11E,dyn_m11F,dyn_m11W);
	dyn_m12A=Setup_GetFirstIntFromLine(fin);
	dyn_m12B=Setup_GetFirstIntFromLine(fin);
	dyn_m12D=Setup_GetFirstIntFromLine(fin);
	dyn_m12E=Setup_GetFirstIntFromLine(fin);
	dyn_m12K=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m12A %d dyn_m12B %d dyn_m12D %d dyn_m12E %d dyn_m12K %d\n",dyn_m12A,dyn_m12B,dyn_m12D,dyn_m12E,dyn_m12K);
	dyn_m13A=Setup_GetFirstIntFromLine(fin);
	dyn_m13B=Setup_GetFirstIntFromLine(fin);
	dyn_m13K=Setup_GetFirstIntFromLine(fin);
	printf("dyn_m13A %d dyn_m13B %d dyn_m13K %d\n",dyn_m13A,dyn_m13B,dyn_m13K);
	dyn_mFCC=Setup_GetFirstIntFromLine(fin);
	dyn_mHCP=Setup_GetFirstIntFromLine(fin);
	dyn_mBCC_9=Setup_GetFirstIntFromLine(fin);
	dyn_mBCC_15=Setup_GetFirstIntFromLine(fin);
	printf("do_mFCC %d do_mHCP %d do_mBCC_9 %d do_mBCC_15 %d\n\n",dyn_mFCC,dyn_mHCP,dyn_mBCC_9,dyn_mBCC_15);

	fclose(fin);
}

void Setup_ReadXmol() { // reads coordinates and particle types from xmol trajectory
	int i, j, k, f, remainder, write;
	char c;
	double tx, ty, tz;
	int cntA;
	char input[1000],errMsg[1000];

	f=0;
	for (i=0; i<FRAMES; i++) {
		remainder=i%SAMPLEFREQ;
		if (remainder==0 && f<FRAMES && i>=STARTFROM) {
			write=1;
		}
		else write=0;
		if (f>=FRAMES) break;

		fscanf(fXmol,"%d\n", &k);
		if (k!=N) {
			sprintf(errMsg,"Setup_ReadXmol(): N %d from input frame %d does not match N %d from params file\n",k,i,N);
			Error(errMsg);
		}
		fgets(input,1000,fXmol);
		cntA=0;
		for (j=0; j<N; ++j) {
			if (feof(fXmol)) Error("Setup_ReadXmol(): end of input file reached\n");
			fscanf(fXmol,"	%c	%lg	%lg	%lg\n", &c,&tx,&ty,&tz);
			if (c=='A') cntA++;
			else if (c=='C') cntA++;
			if (write==1) {
				if (c=='A') rtype[j]=1;
				else if (c=='B') rtype[j]=2;
				else if (c=='C') rtype[j]=1;
				else if (c=='D') rtype[j]=2;
				else {
					sprintf(errMsg,"Setup_ReadXmol(): unrecognized character of particle j %d from input frame %d\n",j,i);
					Error(errMsg);
				}

				r[f][0][j]=tx;	r[f][1][j]=ty;	r[f][2][j]=tz;
				if (PRINTINFO==1) if (j==N-1) printf("f%d part%d %c %.5lg %.5lg %.5lg\n\n",i,j,c,r[f][0][j],r[f][1][j],r[f][2][j]);
			}
		}
		if (cntA!=NA) {
			sprintf(errMsg,"Setup_ReadXmol(): NA %d from input frame %d does not match NA %d from params file\n",cntA,i,NA);
			Error(errMsg);
		}
		if (write==1) f++;
	}

	fclose(fXmol);
}

void Setup_ReadBonds() { // reads bonds in from bonds trajectory file
	int f, i, j, k, remainder;
	int part_index, no_bonds;
	int bonded_part;
	double bond_length;
	char input[1000], errMsg[1000];


	f=0;
	for (i=0; i<FRAMES; i++) {
		remainder=i%SAMPLEFREQ;
		if (remainder==0 && f<FRAMES) {
			fgets(input,1000,fBonds);
			for (j=0; j<N; ++j) {
				fscanf(fBonds,"%d	%d", &part_index, &no_bonds);
				if (j!=part_index) {
					sprintf(errMsg,"Setup_ReadBonds(): part_index %d from input frame %d  does not match j %d\n",part_index,f,j);
					Error(errMsg);
				}
				if (no_bonds>nB) {
					sprintf(errMsg,"Setup_ReadBonds(): no_bonds %d part %d from input frame %d is greater than nB %d\n",no_bonds,part_index,f,nB);
					Error(errMsg);
				}
				cnb[f][j]=no_bonds;
				if (no_bonds>0) {
					for (k=0; k<no_bonds; ++k) {
						fscanf(fBonds,"	%d	%lg", &bonded_part, &bond_length);
						bNums[f][j][k]=bonded_part;
					}
				}
				fscanf(fBonds,"\n");
			}
			f++;
		}
	}

	fclose(fBonds);
}

int Bonds_BondCheck(int f, int i, int j) {	// Returns 1 if i & j are bonded; 0 otherwise
	int k;

	for (k=0; k<cnb[f][i]; ++k) {
		if (bNums[f][i][k] == j) return 1;
	}
	return 0;
}

void Setup_InitVars() { // Initialize lots of important variables
	int f, i, j, k;
	char errMsg[1000];

	nobins=(int)(0.5*sidex/binWidth);
	range2=nobins*binWidth;
	range2=range2*range2;

	fXmol=fopen(fXmolName,"r"); // xmol trajectory open
	if (fXmol==NULL)  {
		sprintf(errMsg,"Setup_InitVars() : Error opening file %s",fXmolName);	// Always test file open
		Error(errMsg);
	}

	fBonds=fopen(fBondsName,"r"); // bonds trajectory open
	if (fBonds==NULL)  {
		sprintf(errMsg,"Setup_InitVars() : Error opening file %s",fBondsName);	// Always test file open
		Error(errMsg);
	}

	// particle type array
	rtype=malloc(N*sizeof(int)); if (rtype==NULL) { sprintf(errMsg,"Setup_InitVars(): rtype[] malloc out of memory\n");	Error(errMsg); }
	for (k=0;k<N;k++) { // test fill array
		rtype[k]=-1;
	}

	// particle position array
	r=malloc(FRAMES*sizeof(double **));
	if (r==NULL) { sprintf(errMsg,"Setup_InitVars(): r[] malloc out of memory\n");	Error(errMsg); }
	for (i=0;i<FRAMES;i++) {
		r[i]=malloc(3*sizeof(double *));
		if (r[i]==NULL) { sprintf(errMsg,"Setup_InitVars(): r[][] malloc out of memory\n");	Error(errMsg); }
		for (j=0;j<3;j++) {
			r[i][j]=malloc(N*sizeof(double));
			if (r[i][j]==NULL) { sprintf(errMsg,"Setup_InitVars(): r[][][] malloc out of memory\n");	Error(errMsg); }
		}
	}
	for (i=0;i<FRAMES;i++) { // test fill array
		for (j=0;j<3;j++) {
			for (k=0;k<N;k++) {
				r[i][j][k]=-1000.0;
			}
		}
	}

	// number of bonds per particle per frame array
	cnb = malloc(FRAMES*sizeof(int *));	if (cnb==NULL) { sprintf(errMsg,"Setup_InitVars(): cnb[] malloc out of memory\n");	Error(errMsg); }
	for (f=0; f<FRAMES; ++f) { cnb[f] = malloc(N*sizeof(int));	if (cnb[f]==NULL) { sprintf(errMsg,"Setup_InitVars(): cnb[][] malloc out of memory\n");	Error(errMsg); } }
	// bonds per particle per frame array
	bNums = malloc(FRAMES*sizeof(int **));	if (bNums==NULL) { sprintf(errMsg,"Setup_InitVars(): bNums[] malloc out of memory\n");	Error(errMsg); }
	for (f=0; f<FRAMES; ++f) { bNums[f] = malloc(N*sizeof(int *));	if (bNums[f]==NULL) { sprintf(errMsg,"Setup_InitVars(): bNums[][] malloc out of memory\n");	Error(errMsg); } }
	for (f=0; f<FRAMES; ++f) {
		for (i=0; i<N; i++) {
			bNums[f][i] = malloc(nB*sizeof(int));	if (bNums[f][i]==NULL) { sprintf(errMsg,"Setup_InitVars(): bNums[][][] malloc out of memory\n");	Error(errMsg); }
		}
	}
	for (f=0; f<FRAMES; ++f) { // test fill array
		for (i=0; i<N; i++) {
			for (j=0; j<nB; j++) {
				bNums[f][i][j]=-1;
			}
			cnb[f][i]=-1;
		}
	}

	// lengths of bonds in a single frame array
	bondlengths = malloc(N*sizeof(double *));	if (bondlengths==NULL) { sprintf(errMsg,"Setup_InitVars(): bondlengths[] malloc out of memory\n");	Error(errMsg); }
	for (j=0; j<N; ++j) { bondlengths[j] = malloc(nB*sizeof(double));	if (bondlengths[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): bondlengths[][] malloc out of memory\n");	Error(errMsg); } }
	for (i=0; i<N; i++) {
		for (j=0; j<nB; j++) {
			bondlengths[i][j]=-1000.0;
		}
	}

	calc_histo = malloc(FRAMES*sizeof(int));	if (calc_histo==NULL) { sprintf(errMsg,"Setup_InitVars(): calc_histo[] malloc out of memory\n");	Error(errMsg); }
	for (f=0; f<FRAMES; ++f) calc_histo[f]=0;

	mClust=0;

	if (doRegions==1 || doSplits==1) {
		used_part = malloc(N*sizeof(int));	if (used_part==NULL) { sprintf(errMsg,"Setup_InitVars(): used_part[] malloc out of memory\n");	Error(errMsg); }
		used_bond = malloc(N*sizeof(int *));	if (used_bond==NULL) { sprintf(errMsg,"Setup_InitVars(): used_bond[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { used_bond[j] = malloc(nB*sizeof(int));	if (used_bond[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): used_bond[][] malloc out of memory\n");	Error(errMsg); } }
		for (i=0; i<N; i++) {
			used_part[i]=-1;
			for (j=0; j<nB; j++) {
				used_bond[i][j]=-1;
			}
		}
	}

	if (doRegions==1) {
		part_move = malloc(N*sizeof(int *));	if (part_move==NULL) { sprintf(errMsg,"Setup_InitVars(): part_move[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { part_move[j] = malloc(3*sizeof(int));	if (part_move[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): part_move[][] malloc out of memory\n");	Error(errMsg); } }
		bond_type_x = malloc(N*sizeof(int *));	if (bond_type_x==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_x[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { bond_type_x[j] = malloc(nB*sizeof(int));	if (bond_type_x[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_x[][] malloc out of memory\n");	Error(errMsg); } }
		bond_type_y = malloc(N*sizeof(int *));	if (bond_type_y==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_y[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { bond_type_y[j] = malloc(nB*sizeof(int));	if (bond_type_y[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_y[][] malloc out of memory\n");	Error(errMsg); } }
		bond_type_z = malloc(N*sizeof(int *));	if (bond_type_z==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_z[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { bond_type_z[j] = malloc(nB*sizeof(int));	if (bond_type_z[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): bond_type_z[][] malloc out of memory\n");	Error(errMsg); } }
		used_part_noPBCs = malloc(N*sizeof(int));	if (used_part_noPBCs==NULL) { sprintf(errMsg,"Setup_InitVars(): used_part_noPBCs[] malloc out of memory\n");	Error(errMsg); }
		used_bond_noPBCs = malloc(N*sizeof(int *));	if (used_bond_noPBCs==NULL) { sprintf(errMsg,"Setup_InitVars(): used_bond_noPBCs[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<N; ++j) { used_bond_noPBCs[j] = malloc(nB*sizeof(int));	if (used_bond_noPBCs[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): used_bond_noPBCs[][] malloc out of memory\n");	Error(errMsg); } }
		no_regions = malloc(FRAMES*sizeof(int));	if (no_regions==NULL) { sprintf(errMsg,"Setup_InitVars(): no_regions[] malloc out of memory\n");	Error(errMsg); }
		regions = malloc(FRAMES*sizeof(int *));	if (regions==NULL) { sprintf(errMsg,"Setup_InitVars(): regions[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<FRAMES; ++j) { regions[j] = malloc(max_regions*sizeof(int));	if (regions[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): regions[][] malloc out of memory\n");	Error(errMsg); } }
		regions_noPBCs = malloc(max_regions*sizeof(int));	if (regions_noPBCs==NULL) { sprintf(errMsg,"Setup_InitVars(): regions_noPBCs[] malloc out of memory\n");	Error(errMsg); }
		regions_Rg = malloc(FRAMES*sizeof(double *));	if (regions_Rg==NULL) { sprintf(errMsg,"Setup_InitVars(): regions_Rg[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<FRAMES; ++j) { regions_Rg[j] = malloc(max_regions*sizeof(double));	if (regions_Rg[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): regions_Rg[][] malloc out of memory\n");	Error(errMsg); } }
		regions_mean_cluster_lifetime = malloc(FRAMES*sizeof(double *));	if (regions_mean_cluster_lifetime==NULL) { sprintf(errMsg,"Setup_InitVars(): regions_mean_cluster_lifetime[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<FRAMES; ++j) { regions_mean_cluster_lifetime[j] = malloc(max_regions*sizeof(double));	if (regions_mean_cluster_lifetime[j]==NULL) { sprintf(errMsg,"Setup_InitVars(): regions_mean_cluster_lifetime[][] malloc out of memory\n");	Error(errMsg); } }

		for (i=0; i<N; i++) {
			part_move[i][0]=part_move[i][1]=part_move[i][2]=0;
			for (j=0; j<nB; j++) {
				bond_type_x[i][j]=bond_type_y[i][j]=bond_type_z[i][j]=0;
			}
			used_part_noPBCs[i]=-1;
			for (j=0; j<nB; j++) {
				used_bond_noPBCs[i][j]=-1;
			}
		}
		for (i=0; i<FRAMES; i++) {
			for (j=0; j<max_regions; j++) {
				regions[i][j]=0;
				regions_noPBCs[j]=0;
				regions_Rg[i][j]=0.0;
				regions_mean_cluster_lifetime[i][j]=0.0;
			}
			no_regions[i]=0;
		}
	}

	/*if (output_neighbours==1) {
		temp_neighbours = malloc(FRAMES*sizeof(int *));	if (temp_neighbours==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_neighbours[] malloc out of memory\n");	Error(errMsg); }
		for (f=0; f<FRAMES; ++f) {
			temp_neighbours[f] = malloc(clust_nB*sizeof(int));
			if (temp_neighbours[f]==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_neighbours[][] malloc out of memory\n");
				Error(errMsg);
			}
		}

		no_temp_neighbours = malloc(FRAMES*sizeof(int));	if (no_temp_neighbours==NULL) { sprintf(errMsg,"Input_Cluster_Type(): no_temp_neighbours[] malloc out of memory\n");	Error(errMsg); }

		for (f=0; f<FRAMES; ++f) {
			for (j=0; j<clust_nB; j++) {
				temp_neighbours[f][j]=-1;
			}
			no_temp_neighbours[f]=0;
		}
	}*/
}

void Setup_ResetVars() {	// Reset variables for each frame
	int i, j;

	for (i=0; i<N; ++i) { // reset bond lengths
		for (j=0;j<nB;++j) {
			bondlengths[i][j]=-1000.0;
		}
	}
	if (doRegions==1 || doSplits==1) {
		for (i=0; i<N; ++i) {
			used_part[i]=-1;
			for (j=0; j<nB; j++) {
				used_bond[i][j]=-1;
			}
		}
	}

	if (doRegions==1) { // reset arrays relating to cluster domains analysis
		for (i=0; i<N; i++) {
			part_move[i][0]=part_move[i][1]=part_move[i][2]=0;
			for (j=0; j<nB; j++) {
				bond_type_x[i][j]=bond_type_y[i][j]=bond_type_z[i][j]=0;
			}
			used_part_noPBCs[i]=-1;
			for (j=0; j<nB; j++) {
				used_bond_noPBCs[i][j]=-1;
			}
		}
	}
}

void Setup_FreeVars() {	// Free bond detection variables
	int f, i;

	free(rtype);

	for (f=0; f<FRAMES; ++f) {
		for (i=0; i<3; ++i) free(r[f][i]);
		for (i=0; i<N; ++i) {
			free(bNums[f][i]);
		}
	}
	for (f=0; f<FRAMES; ++f) {
		free(r[f]);
		free(bNums[f]);
		free(cnb[f]);
	}
	free(r);
	free(bNums);
	free(cnb);

	for (i=0; i<N; ++i) {
		free(bondlengths[i]);
	}
	free(bondlengths);

	free(calc_histo);

	if (doRegions==1 || doSplits==1) {
		for (i=0; i<N; ++i) {
			free(used_bond[i]);
		}
		free(used_bond);
		free(used_part);
	}
	if (doRegions==1) {
		for (i=0; i<N; ++i) {
			free(part_move[i]);
			free(bond_type_x[i]);
			free(bond_type_y[i]);
			free(bond_type_z[i]);
		}
		free(bond_type_x);
		free(bond_type_y);
		free(bond_type_z);
		free(part_move);
		for (i=0; i<N; ++i) {
			free(used_bond_noPBCs[i]);
		}
		free(used_bond_noPBCs);
		free(used_part_noPBCs);
		for (i=0; i<FRAMES; ++i) {
			free(regions[i]);
			free(regions_Rg[i]);
			free(regions_mean_cluster_lifetime[i]);
		}
		free(regions);
		free(regions_noPBCs);
		free(regions_Rg);
		free(regions_mean_cluster_lifetime);
		free(no_regions);
	}

	/*if (output_neighbours==1) {
		for (i=0; i<FRAMES; i++) {
			free(temp_neighbours[i]);
		}
		free(temp_neighbours);
		free(no_temp_neighbours);
	}*/
}

int *resize_1D_int(int *the_array, int old_col_size, int new_col_size) {
	int i;
	char errMsg[1000];

	the_array=realloc(the_array,new_col_size*sizeof(int));
	if (the_array == NULL) { sprintf(errMsg,"resize_1D_int(): the_array[] out of memory old_col_size %d new_col_size %d\n",old_col_size,new_col_size); Error(errMsg); }
	for (i=old_col_size; i<new_col_size; i++) the_array[i]=-1;

	return the_array;
}

void Input_Cluster_Type(char *clustName, int clustSize, int subClust, char *strSubClust, int centre_shell) { // read in clusters from dyn_*** file outputted by dynamic tcc
	int throwaway;
	int i, j, sum;
	int maxNoLifetimes;
	char errMsg[1000], input[1000], output[1000], cthrowaway[50000];
	FILE *fIn, *fOut;
	double longest_single_instance, longest_total_instance, longest_max_instance;
	int longest_single_int_instance, longest_total_int_instance, longest_max_int_instance;

	// open dyn_** file
	sprintf(input,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.dyn_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs,clustName);
	fIn=fopen(input, "r");
	if (fIn==NULL)  { // test opened ok
		sprintf(errMsg,"Input_Cluster_Type(): Error opening file %s",input);	// Always test file open
		Error(errMsg);
	}

	fgets(cthrowaway,10000,fIn); // throwaway first information line
	mClust=0; // read in number of clusters that the file contains
	fscanf(fIn,"%d\n",&mClust);
	printf("number of %s clusters coming in %d\n",clustName,mClust);
	if (mClust==0) { // no clusters present
		printf("no clusters in dynamic tcc file %s\n",input);
		fclose(fIn);
		return;
	}
	fgets(cthrowaway,50000,fIn); // throw away column headings line

	// set up array for particles in each unique cluster
	temp_parts = malloc(mClust*sizeof(int *));	if (temp_parts==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_parts[] malloc out of memory\n");	Error(errMsg); }
	for (j=0; j<mClust; ++j) { temp_parts[j] = malloc(clustSize*sizeof(int));	if (temp_parts[j]==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_parts[][] malloc out of memory\n");	Error(errMsg); } }

	// set up array for number of A type particles in each unique cluster
	if (doBinary==1) {
		temp_parts_A = malloc(mClust*sizeof(int));	if (temp_parts_A==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_parts_A[] malloc out of memory\n");	Error(errMsg); }
		if (centre_shell==1) { // set up array for number of A type particles as centre or in shell of each unique cluster
			temp_parts_A_cen = malloc(mClust*sizeof(int));	if (temp_parts_A_cen==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_parts_A_cen[] malloc out of memory\n");	Error(errMsg); }
			temp_parts_A_shell = malloc(mClust*sizeof(int));	if (temp_parts_A_shell==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_parts_A_shell[] malloc out of memory\n");	Error(errMsg); }
		}
	}

	// set up array for tcc lifetimes of each unique cluster
	temp_lifetimes = malloc(mClust*sizeof(double *));	if (temp_lifetimes==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_lifetimes[] malloc out of memory\n");	Error(errMsg); }
	for (j=0; j<mClust; ++j) { temp_lifetimes[j] = malloc(4*sizeof(double));	if (temp_lifetimes[j]==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_lifetimes[][] malloc out of memory\n");	Error(errMsg); } }

	// set up array for number of, starting frames, ending frames for events of each unique cluster ### POOR PERFORMANCE COS VERY BIG ARRAYS
	temp_events = malloc(mClust*sizeof(int *));	if (temp_events==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_events[] malloc out of memory\n");	Error(errMsg); }
	for (j=0; j<mClust; ++j) { temp_events[j] = malloc(1*sizeof(int));	if (temp_events[j]==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_events[][] malloc out of memory\n");	Error(errMsg); } }

	if (doSubClusts==1) { // set up arrays for the unique id's of subclusters, if they exist
		temp_sub = malloc(mClust*sizeof(int *));	if (temp_sub==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_sub[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<mClust; ++j) { temp_sub[j] = malloc(subClust*sizeof(int));	if (temp_sub[j]==NULL) { sprintf(errMsg,"Input_Cluster_Type(): temp_sub[][] malloc out of memory\n");	Error(errMsg); } }
	}

	if (doSplits==1) { // make sure no subcluster of a cluster becomes unbonded during a disappearence event less than or equal to go_away_for
		temp_splits = malloc(mClust*sizeof(int *));	if (temp_splits==NULL) { sprintf(errMsg,"main(): temp_splits[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<mClust; ++j) { temp_splits[j] = malloc(1*sizeof(int));	if (temp_splits[j]==NULL) { sprintf(errMsg,"main(): temp_splits[][] malloc out of memory\n");	Error(errMsg); } }
		for (i=0; i<mClust; ++i) {
			for (j=0; j<1; ++j) {
				temp_splits[i][j]=-1;
			}
		}
	}

	for (i=0; i<mClust; i++) { // test write to all of the arrays
		for (j=0; j<clustSize; j++) {
			temp_parts[i][j]=-1;
		}
		if (doBinary==1) {
			temp_parts_A[i]=-1;
			if (centre_shell==1) {
				temp_parts_A_cen[i]=-1;
				temp_parts_A_shell[i]=-1;
			}
		}
		for (j=0; j<4; j++) {
			temp_lifetimes[i][j]=-1.0;
		}
		for (j=0; j<1; j++) {
			temp_events[i][j]=-1;
		}
		if (doSubClusts==1) {
			for (j=0; j<subClust; j++) {
				temp_sub[i][j]=-1;
			}
		}
	}

	maxNoLifetimes=0;
	longest_single_instance=longest_total_instance=longest_max_instance=0.0;
	longest_single_int_instance=longest_total_int_instance=longest_max_int_instance=0; // find longest living cluster
	for (i=0; i<mClust; i++) { // loop over all input clusters
		fscanf(fIn,"%d	%lg	%lg	%lg	%lg",&throwaway,&temp_lifetimes[i][0],&temp_lifetimes[i][1],&temp_lifetimes[i][2],&temp_lifetimes[i][3]);
		if (longest_single_instance<temp_lifetimes[i][0]) longest_single_instance=temp_lifetimes[i][0];
		if (longest_total_instance<temp_lifetimes[i][1]) longest_total_instance=temp_lifetimes[i][1];
		if (longest_max_instance<temp_lifetimes[i][2]) longest_max_instance=temp_lifetimes[i][2];
		if (doBinary==1) {
			fscanf(fIn,"	%d",&temp_parts_A[i]);
			if (centre_shell==1) fscanf(fIn,"	%d	%d",&temp_parts_A_cen[i],&temp_parts_A_shell[i]);
		}
		for (j=0; j<clustSize; j++) {
			fscanf(fIn,"	%d",&temp_parts[i][j]);
		}
		if (doSubClusts==1) {
			for (j=0; j<subClust; j++) fscanf(fIn,"	%d",&temp_sub[i][j]);
		}
		fscanf(fIn,"	%d",&temp_events[i][0]);
		if (maxNoLifetimes<temp_events[i][0]) maxNoLifetimes=temp_events[i][0];
		temp_events[i]=resize_1D_int(temp_events[i], 1, 1+2*temp_events[i][0]);
		if (doSplits==1) temp_splits[i]=resize_1D_int(temp_splits[i], 1, 1+2*temp_events[i][0]);
		sum=0;
		for (j=0; j<temp_events[i][0]; j++) {
			fscanf(fIn,"	%d	%d",&temp_events[i][2*j+1],&temp_events[i][2*j+2]);
			if (longest_single_int_instance<temp_events[i][2*j+2]-temp_events[i][2*j+1]) longest_single_int_instance=temp_events[i][2*j+2]-temp_events[i][2*j+1];
			sum+=temp_events[i][2*j+2]-temp_events[i][2*j+1];
		}
		if (longest_total_int_instance<sum) longest_total_int_instance=sum;
		if (longest_max_int_instance<temp_events[i][2*temp_events[i][0]]-temp_events[i][1]) longest_max_int_instance=temp_events[i][2*temp_events[i][0]]-temp_events[i][1];
		fscanf(fIn,"\n");
	}
	fclose(fIn);
	printf("\nRead %s\n",input);
	printf("by frame count longest single total max lifetime\n");
	printf("%d %d %d\n",longest_single_int_instance,longest_total_int_instance,longest_max_int_instance);

	printf("by t_alpha longest single total max lifetimes\n");
	printf("%lg %lg %lg\n",longest_single_instance,longest_total_instance,longest_max_instance);

	if (doDynFileRewrite==1) { // rewrite dyn_** file (just for debugging - test to see it was all read in ok)
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.test.dyn_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs,clustName);
		fOut=fopen(output, "w");
		if (fOut==NULL)  {
			sprintf(errMsg,"Input_Cluster_Type() : Error opening file %s",output);	// Always test file open
			Error(errMsg);
		}

		fprintf(fOut,"%s\n",output);
		fprintf(fOut,"%d\n",mClust);
		fprintf(fOut,"i	longest_single_lifetime	total_lifetime	max_lifetime	ratio");
		if (doBinary==1) {
			fprintf(fOut,"	N_A");
			if (centre_shell==1) fprintf(fOut,"	N_A_cen	N_A_shell");
		}
		for (j=0;j<clustSize;j++) fprintf(fOut,"	p%d",j+1);
		if (doSubClusts==1) fprintf(fOut,"	%s",strSubClust);
		fprintf(fOut,"	events");
		for (j=0;j<maxNoLifetimes;j++) fprintf(fOut,"	l_s_%d	l_e_%d",j,j);
		fprintf(fOut,"\n");
		for (i=0;i<mClust;i++) {
			fprintf(fOut,"%d	%.15lg	%.15lg	%.15lg	%.15lg",i,temp_lifetimes[i][0],temp_lifetimes[i][1],temp_lifetimes[i][2],temp_lifetimes[i][3]);
			if (doBinary==1) {
				fprintf(fOut,"	%d",temp_parts_A[i]);
				if (centre_shell==1) fprintf(fOut,"	%d	%d",temp_parts_A_cen[i],temp_parts_A_shell[i]);
			}
			for (j=0;j<clustSize;j++) fprintf(fOut,"	%d",temp_parts[i][j]);
			if (doSubClusts==1) {
				for (j=0;j<subClust;j++) fprintf(fOut,"	%d",temp_sub[i][j]);
			}
			fprintf(fOut,"	%d",temp_events[i][0]);
			for (j=1;j<2*temp_events[i][0]+1;j++) fprintf(fOut,"	%d",temp_events[i][j]);
			fprintf(fOut,"\n");
		}
		fclose(fOut);
		printf("Rewritten dyn_... file to check all is ok %s\n",output);
	}
}

void Free_temp(int centre_shell) { // free the temp variables after the new clusters have been found
	int i;
	for (i=0; i<mClust; i++) {
		free(temp_parts[i]);
		free(temp_lifetimes[i]);
		free(temp_events[i]);
		if (doSubClusts==1) free(temp_sub[i]);
	}
	free(temp_parts);
	free(temp_lifetimes);
	free(temp_events);
	if (doSubClusts==1) free(temp_sub);

	if (doBinary==1) {
		free(temp_parts_A);
		if (centre_shell==1) {
			free(temp_parts_A_cen);
			free(temp_parts_A_shell);
		}
	}
}

void DFS_loop2(int f, int in_part) {
	int j, k;

	used_part[in_part]=1;
	for (j=0; j<cnb[f][in_part]; j++) {
		if (used_bond[in_part][j]!=0) continue;
		if (used_part[bNums[f][in_part][j]]!=0) continue;
		used_bond[in_part][j]=1;
		for (k=0; k<cnb[f][bNums[f][in_part][j]]; k++) {
			if (bNums[f][bNums[f][in_part][j]][k]!=in_part) continue;
			used_bond[bNums[f][in_part][j]][k]=1;
			break;
		}
		DFS_loop2(f,bNums[f][in_part][j]);
	}
}

int No_Break_Up(int clustSize, int *clust_parts, int f) {
	int i, j, k;
	int regioncnt;

	for (i=0; i<N; i++) {
		used_part[i]=-1;
		for (j=0; j<nB; j++) {
			used_bond[i][j]=-1;
		}
	}

	for (i=0; i<clustSize; i++) {
		used_part[clust_parts[i]]=0;
		for (j=0; j<cnb[f][clust_parts[i]]; j++) {
			for (k=0; k<clustSize; k++) {
				if (clust_parts[i]==clust_parts[k]) continue;
				if (bNums[f][clust_parts[i]][j]==clust_parts[k]) used_bond[clust_parts[i]][j]=0;
			}
		}
	}

	regioncnt=0;
	for (i=0; i<clustSize; i++) {
		if (used_part[clust_parts[i]]!=0) continue;
		DFS_loop2(f,clust_parts[i]);
		regioncnt++;
	}

	return regioncnt;
}

void Find_Splits(int clustSize, int no_of_clusters) {
	int e, f, i, j, k;
	int *clust_parts;
	int break_up, between_frame, remainder;
	int write;
	char errMsg[1000];

	clust_parts = malloc(clustSize*sizeof(int));	if (clust_parts==NULL) { sprintf(errMsg,"Find_Splits(): clust_parts[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<clustSize; i++) clust_parts[i]=-1;

	f=0;
	for (e=0;e<FRAMES;e++) {
		if (f==FRAMES) break;
		remainder=e%SAMPLEFREQ;
		if (remainder==0 && f<FRAMES && e>=STARTFROM) {
			write=1;
		}
		else write=0;
		if (write==1) {
			for (i=0; i<no_of_clusters; i++) {
				between_frame=0;
				if (temp_events[i][0]==1) continue; // no splits of this cluster as no disappearence events
				for (j=0; j<temp_events[i][0]-1; j++) {	// loop over all events to look at disappearences
					if (temp_splits[i][j]==-1) temp_splits[i][j]=0;
					if (temp_events[i][2*j+3]-temp_events[i][2*j+2]>go_away_for) {
						between_frame=0;
						continue;
					}
					if (temp_splits[i][j]==1) { // already seen a split before in a disappearnce event of this cluster of this cluster
						between_frame=0;
						continue;
					}
					if (f>temp_events[i][2*j+2] && f<temp_events[i][2*j+3]) {
						between_frame=1;
						break;
					}
				}
				if (between_frame==0) continue;
				for (k=0; k<clustSize; k++) {
					clust_parts[k]=temp_parts[i][k];
				}
				break_up=No_Break_Up(clustSize,&clust_parts[0],f);
				if (break_up>1) temp_splits[i][j]=1;
			}
			f++;
		}
	}

	free(clust_parts);
}

void Alloc_the_clusts(int calc_nosamples, int clustSize) {
	int i,j;
	char errMsg[1000];

	the_clusts = malloc(calc_nosamples*sizeof(int *));	if (the_clusts==NULL) { sprintf(errMsg,"Alloc_the_clusts(): the_clusts[] malloc out of memory\n");	Error(errMsg); }
	for (j=0; j<calc_nosamples; ++j) { the_clusts[j] = malloc((clustSize+3)*sizeof(int));	if (the_clusts[j]==NULL) { sprintf(errMsg,"Alloc_the_clusts(): the_clusts[][] malloc out of memory\n");	Error(errMsg); } }
	if (doBinary==1) {
		the_clusts_mixture = malloc(calc_nosamples*sizeof(int *));	if (the_clusts_mixture==NULL) { sprintf(errMsg,"Alloc_the_clusts(): the_clusts_mixture[] malloc out of memory\n");	Error(errMsg); }
		for (j=0; j<calc_nosamples; ++j) { the_clusts_mixture[j] = malloc(3*sizeof(int));	if (the_clusts_mixture[j]==NULL) { sprintf(errMsg,"Alloc_the_clusts(): the_clusts_mixture[][] malloc out of memory\n");	Error(errMsg); } }
	}
	for (i=0; i<calc_nosamples; ++i) {
		for (j=0; j<clustSize+3; ++j) {
			the_clusts[i][j]=-1;
		}
		if (doBinary==1) {
			for (j=0; j<3; ++j) {
				the_clusts_mixture[i][j]=-1;
			}
		}
	}
	clustered_parts = malloc(FRAMES*sizeof(int *));	if (clustered_parts==NULL) { sprintf(errMsg,"Alloc_the_clusts(): clustered_parts[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<FRAMES; i++) { clustered_parts[i] = malloc(N*sizeof(int *));	if (clustered_parts[i]==NULL) { sprintf(errMsg,"Alloc_the_clusts(): clustered_parts[][] malloc out of memory\n");	Error(errMsg); } }
	for (i=0; i<FRAMES; i++) {
		for (j=0; j<N; j++) clustered_parts[i][j]=0;
	}
}

void Free_the_clusts(int calc_nosamples) {
	int i;
	for (i=0; i<calc_nosamples; ++i) {
		free(the_clusts[i]);
		if (doBinary==1) free(the_clusts_mixture[i]);
	}
	free(the_clusts);
	if (doBinary==1) free(the_clusts_mixture);

	for (i=0; i<FRAMES; ++i) {
		free(clustered_parts[i]);
	}
	free(clustered_parts);
}

int Calc_Lifetimes(char *fileName, int n_clusters, int clustSize, int binary_center) {
	/* The aim of this function is to post process the clusters outputted by the dynamic TCC.
	A couple of things needs doing to the dynamic TCC clusters to get sensible results:

		1. Remove clusters from first go_away_for and last go_away_for frames of the trajectory.
		   This is because if a clusters is first or last seen in the start or end of the trajectory
		   we can't be sure that it doesnt exist at earlier or later times.
		2. Remove disappearences of a cluster if the disappearence is less than or equal
		    to go_away_for frames and no subset of the particles become unbonded over the disappearence.
	*/

	FILE *file;
	int i, j, k;
	int running, alive_for, break_up;
	int curr_start, curr_end, dropped, calc_nosamples, calc_nosamples_again;
	char errMsg[1000];

	dropped=alive_for=0;
	calc_maxlength=calc_nosamples=0;

	for (i=0; i<FRAMES; i++) {
		calc_histo[i]=0;
	}

	if (doDynamicAnalysis==1) { // remove disappearences with lengths less than or equal to go_away_for
	for (i=0; i<n_clusters; i++) { // loop over all the clusters from dynamic TCC dyn_** file
		running=0;
		curr_start=curr_end=-1;
		for (k=0; k<temp_events[i][0]; k++) {	// loop over all events of ith cluster
			if (temp_events[i][2*k+1]<start_from || temp_events[i][2*k+2]>end_by) { // if event of ith cluster starts (or ends) in first (or last) go_away_for frames drop cluster completely
				dropped++;
				if (running==0) continue;
				else { // running and curr_end event is in end go_away_for frames for trajectory, therefore drop cluster
					running=0;
					break;
				}
			}
			if (running==0) { // set up new lifetime event for ith cluster
				curr_start=temp_events[i][2*k+1];
				curr_end=temp_events[i][2*k+2];
				running=1;
				if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
					calc_histo[curr_end-curr_start]++;
					calc_nosamples++;
					if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
					running=0;
				}
			}
			else { // already running an event of ith cluster, see it it continues over next disappearence
				if (temp_events[i][2*k+1]-curr_end>go_away_for) { // breaks up over disappearence
					calc_histo[curr_end-curr_start]++; // add event to histogram
					calc_nosamples++;
					if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
					running=0;

					curr_start=temp_events[i][2*k+1]; // start new event for ith cluster
					curr_end=temp_events[i][2*k+2];
					running=1;
					if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
						calc_histo[curr_end-curr_start]++;
						calc_nosamples++;
						if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
						running=0;
					}
				}
				else { // disappearence is for less than or equal to go_away_for frames, need to check no subcluster becomes unbonded during this disappearence event
					if (doSplits==1) {
						if (temp_splits[i][k-1]==1) break_up=2;  // already seen split of cluster over this disappearence event
						else break_up=1;
						if (break_up==1) { // no split of cluster over break up event, shuffle end along to end of next event
							curr_end=temp_events[i][2*k+2];
							if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
								calc_histo[curr_end-curr_start]++;
								calc_nosamples++;
								if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
								running=0;
							}
						}
						else { // split of cluster over disappearence event go_away_frames or less
							calc_histo[curr_end-curr_start]++;
							calc_nosamples++;
							if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
							running=0;

							curr_start=temp_events[i][2*k+1];
							curr_end=temp_events[i][2*k+2];
							running=1;

							if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
								calc_histo[curr_end-curr_start]++;
								calc_nosamples++;
								if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
								running=0;
							}
						}
					}
					else { // cluster continues over the disappearcence event of go_away_for or fewer frames
						curr_end=temp_events[i][2*k+2];
						if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
							calc_histo[curr_end-curr_start]++;
							calc_nosamples++;
							if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
							running=0;
						}
					}
				}
			}
		}
	}
	printf("number of clusts in %d no clusts out %d dropped %d useable_frames %d \n",n_clusters,calc_nosamples,dropped,useable_frames);

	Alloc_the_clusts(calc_nosamples,clustSize); // alloc memory for the new dynamic clusters
	printf("allocated memory for the processed dynamic clusters\n");

	calc_nosamples_again=0;
	for (i=0; i<n_clusters; i++) { // need to run above again to put new_clusters into arrays
		running=0;
		curr_start=curr_end=-1;
		for (k=0; k<temp_events[i][0]; k++) { // loop over all events of ith cluster
			if (temp_events[i][2*k+1]<start_from || temp_events[i][2*k+2]>end_by) { // if event of ith cluster starts (or ends) in first (or last) go_away_for frames drop cluster completely
				dropped++;
				if (running==0) continue;
				else { // running and curr_end event is in end go_away_for frames for trajectory, therefore drop cluster
					running=0;
					break;
				}
			}
			if (running==0) {  // set up new lifetime event for ith cluster
				curr_start=temp_events[i][2*k+1];
				curr_end=temp_events[i][2*k+2];
				alive_for=curr_end-curr_start+1;
				running=1;
				if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
					the_clusts[calc_nosamples_again][0]=curr_start;
					the_clusts[calc_nosamples_again][1]=curr_end;
					if (doBinary==1) {
						the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
						if (binary_center==1) {
							the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
							the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
						}
					}
					for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
					the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
					calc_nosamples_again++;
					running=0;
				}
			}
			else {  // already running an event of ith cluster, see it it continues over next disappearence
				if (temp_events[i][2*k+1]-curr_end>go_away_for) {
					the_clusts[calc_nosamples_again][0]=curr_start;
					the_clusts[calc_nosamples_again][1]=curr_end;
					if (doBinary==1) {
						the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
						if (binary_center==1) {
							the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
							the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
						}
					}
					for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
					the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
					calc_nosamples_again++;
					running=0;

					curr_start=temp_events[i][2*k+1];
					curr_end=temp_events[i][2*k+2];
					alive_for=curr_end-curr_start+1;
					running=1;
					if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
						the_clusts[calc_nosamples_again][0]=curr_start;
						the_clusts[calc_nosamples_again][1]=curr_end;
						if (doBinary==1) {
							the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
							if (binary_center==1) {
								the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
								the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
							}
						}
						for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
						the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
						calc_nosamples_again++;
						running=0;
					}
				}
				else { // disappearence is for less than or equal to go_away_for frames, need to check no subcluster becomes unbonded during this disappearence event
					if (doSplits==1) {
						if (temp_splits[i][k-1]==1) break_up=2;  // already seen split of cluster over this disappearence event
						else break_up=1;
						if (break_up==1) { // no split of cluster over break up event, shuffle end along to end of next event
							curr_end=temp_events[i][2*k+2];
							alive_for+=curr_end-temp_events[i][2*k+1]+1;
							if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
								the_clusts[calc_nosamples_again][0]=curr_start;
								the_clusts[calc_nosamples_again][1]=curr_end;
								if (doBinary==1) {
									the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
									if (binary_center==1) {
										the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
										the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
									}
								}
								for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
								the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
								calc_nosamples_again++;
								running=0;
							}
						}
						else { // split of cluster over disappearence event go_away_frames or less
							the_clusts[calc_nosamples_again][0]=curr_start;
							the_clusts[calc_nosamples_again][1]=curr_end;
							if (doBinary==1) {
								the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
								if (binary_center==1) {
									the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
									the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
								}
							}
							for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
							the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
							calc_nosamples_again++;
							running=0;

							curr_start=temp_events[i][2*k+1];
							curr_end=temp_events[i][2*k+2];
							alive_for=curr_end-curr_start+1;
							running=1;
							if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
								the_clusts[calc_nosamples_again][0]=curr_start;
								the_clusts[calc_nosamples_again][1]=curr_end;
								if (doBinary==1) {
									the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
									if (binary_center==1) {
										the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
										the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
									}
								}
								for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
								the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
								calc_nosamples_again++;
								running=0;
							}
						}
					}
					else { // cluster continues over the disappearcence event of go_away_for or fewer frames
						curr_end=temp_events[i][2*k+2];
						alive_for+=curr_end-temp_events[i][2*k+1]+1;
						if (k==temp_events[i][0]-1) { // if its the last event of ith cluster add it to the histogram
							the_clusts[calc_nosamples_again][0]=curr_start;
							the_clusts[calc_nosamples_again][1]=curr_end;
							if (doBinary==1) {
								the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
								if (binary_center==1) {
									the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
									the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
								}
							}
							for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
							the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
							calc_nosamples_again++;
							running=0;
						}
					}
				}
			}
		}
	}
	if (calc_nosamples_again!=calc_nosamples) Error("calc_nosamples_again!=calc_nosamples");
	}

	else {
	for (i=0; i<n_clusters; i++) {
		running=0;
		curr_start=curr_end=-1;
		for (k=0; k<temp_events[i][0]; k++) {
			if (temp_events[i][2*k+1]<start_from || temp_events[i][2*k+2]>end_by) {
				dropped++;
				continue;
			}

			curr_start=temp_events[i][2*k+1];
			curr_end=temp_events[i][2*k+2];

			calc_histo[curr_end-curr_start]++;
			calc_nosamples++;
			if (calc_maxlength<curr_end-curr_start) calc_maxlength=curr_end-curr_start;
		}
	}
	printf("number of clusts in %d no clusts out %d dropped %d useable_frames %d \n",n_clusters,calc_nosamples,dropped,useable_frames);

	Alloc_the_clusts(calc_nosamples,clustSize); // alloc memory for the new dynamic clusters

	calc_nosamples_again=0;
	for (i=0; i<n_clusters; i++) { // need to run above again to put new_clusters into arrays
		running=0;
		curr_start=curr_end=-1;
		for (k=0; k<temp_events[i][0]; k++) {
			if (temp_events[i][2*k+1]<start_from || temp_events[i][2*k+2]>end_by) {
				dropped++;
				continue;
			}
			curr_start=temp_events[i][2*k+1];
			curr_end=temp_events[i][2*k+2];
			alive_for=curr_end-curr_start+1;
			the_clusts[calc_nosamples_again][0]=curr_start;
			the_clusts[calc_nosamples_again][1]=curr_end;
			if (doBinary==1) {
				the_clusts_mixture[calc_nosamples_again][0]=temp_parts_A[i];
				if (binary_center==1) {
					the_clusts_mixture[calc_nosamples_again][1]=temp_parts_A_cen[i];
					the_clusts_mixture[calc_nosamples_again][2]=temp_parts_A_shell[i];
				}
			}
			for (j=0; j<clustSize; j++) the_clusts[calc_nosamples_again][j+2]=temp_parts[i][j];
			the_clusts[calc_nosamples_again][clustSize+2]=alive_for;
			calc_nosamples_again++;
		}
	}
	if (calc_nosamples_again!=calc_nosamples) Error("calc_nosamples_again!=calc_nosamples");
	}

	printf("filled processed dynamic cluster arrays\n");

	file=fopen(fileName, "w");
	if (file==NULL)  {
		sprintf(errMsg,"Calc_Lifetimes(): Error opening file %s",fileName);	// Always test file open
		Error(errMsg);
	}

	if (doWriteNewClustsFile==1) {
		fprintf(file,"%s\n",fileName);
		fprintf(file,"i	start	end	NA");
		if (binary_center==1) fprintf(file,"	NA cen	NA shell");
		for (i=0; i<clustSize; i++) {
			fprintf(file,"	p%d",i+1);
		}
		fprintf(file,"	alive_for	total length	ratio\n");
		for (i=0; i<calc_nosamples; i++) {
			fprintf(file,"%d	%d	%d",i,the_clusts[i][0],the_clusts[i][1]);
			if (doBinary==1) {
				fprintf(file,"	%d",the_clusts_mixture[i][0]);
				if (binary_center==1) {
					fprintf(file,"	%d",the_clusts_mixture[i][1]);
					fprintf(file,"	%d",the_clusts_mixture[i][2]);
				}
			}
			for (j=0; j<clustSize; j++) 	fprintf(file,"	%d",the_clusts[i][j+2]);
			fprintf(file,"	%d	%d	%.15lg\n",the_clusts[i][clustSize+2],the_clusts[i][1]-the_clusts[i][0]+1,(double)the_clusts[i][clustSize+2]/(the_clusts[i][1]-the_clusts[i][0]+1));
		}

		fclose(file);
	}

	return calc_nosamples;
}

void Write_Lifetime_Distro(char *fileName, int clustSize, int calc_nosamples, int use_binary, int this_NA, int this_NA_cen, int this_NA_shell, int* **in) {
	// writes cluster lifetime correlation functions and Nclus(tau_l)/N
	FILE *file;
	char errMsg[1000];
	int *used_clus;
	int i, j, k, l;
	int *total_histo, *percent, **in_2;
	int calc_nosamples2;
	double *norm_histo, *norm_percent;

	in_2 = malloc(FRAMES*sizeof(int *));	if (in_2==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): in_2[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<FRAMES; i++) { in_2[i] = malloc(N*sizeof(int));	if (in_2[i]==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): in_2[][] malloc out of memory\n");	Error(errMsg); } }
	for (i=0; i<FRAMES; i++) {
		for (j=0; j<N; j++) in_2[i][j]=0;
	}

	used_clus = malloc(calc_nosamples*sizeof(int));	if (used_clus==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): used_clus[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<calc_nosamples; i++) used_clus[i]=0;
	total_histo = malloc(useable_frames*sizeof(int));	if (total_histo==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): total_histo[] malloc out of memory\n");	Error(errMsg); }
	norm_histo = malloc(useable_frames*sizeof(double));	if (norm_histo==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): norm_histo[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<useable_frames; i++) {
		total_histo[i]=0;
		norm_histo[i]=0.0;
	}

	percent = malloc(useable_percent*sizeof(int));	if (percent==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): percent[] malloc out of memory\n");	Error(errMsg); }
	norm_percent = malloc(useable_percent*sizeof(double));	if (norm_percent==NULL) { sprintf(errMsg,"Write_Lifetime_Distro(): norm_percent[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<useable_percent; i++) {
		percent[i]=0;
		norm_percent[i]=0.0;
	}

	calc_nosamples2=0;
	for (i=0; i<FRAMES; i++) {
		calc_histo[i]=0;
	}


	for (j=0; j<calc_nosamples; j++) {
		if (doBinary==1 && use_binary>0) {
			if (the_clusts_mixture[j][0]!=this_NA) continue;
			if (use_binary>1) {
				if (the_clusts_mixture[j][1]!=this_NA_cen) continue;
				if (the_clusts_mixture[j][2]!=this_NA_shell) continue;
			}
		}
		calc_histo[the_clusts[j][1]-the_clusts[j][0]]++;
		calc_nosamples2++;
	}
	total_histo[useable_frames-1]=calc_histo[useable_frames-1];
	for (i=useable_frames-2; i>=0; i--) {
		total_histo[i]+=total_histo[i+1]+calc_histo[i];
	}

	for (i=0; i<useable_frames; i++) {
		norm_histo[i]=(double)total_histo[i]/(double)total_histo[0];
	}

	calc_nosamples2=0;
	for (i=useable_frames-1; i>=0; i--) {
		for (j=0; j<calc_nosamples; j++) {
			if (used_clus[j]==1) continue;
			if (doBinary==1 && use_binary>0) {
				if (the_clusts_mixture[j][0]!=this_NA) { used_clus[j]=1; continue; }
				if (use_binary>1) {
					if (the_clusts_mixture[j][1]!=this_NA_cen) { used_clus[j]=1; continue; }
					if (the_clusts_mixture[j][2]!=this_NA_shell) { used_clus[j]=1; continue; }
				}
			}
			if ((the_clusts[j][1]-the_clusts[j][0])==i) {
				for (k=the_clusts[j][0]; k<=the_clusts[j][1]; k++) {
					for (l=0; l<clustSize; l++) {
						(*in)[k][the_clusts[j][l+2]]=1;
						if (in_2[k][the_clusts[j][l+2]]<i) in_2[k][the_clusts[j][l+2]]=i;
					}
				}
				used_clus[j]=1;
				calc_nosamples2++;
			}
		}
	}

	for (i=useable_percent-1; i>=0; i--) {
		for (j=percent_start; j<=percent_end; j++) {
			for (k=0; k<N; k++) {
				if ((*in)[j][k]==1 && in_2[j][k]>=i) {
					percent[i]+=1;
				}
			}
		}
	}

	for (i=0; i<useable_percent; i++) {
		norm_percent[i]=(double)percent[i]/(N*useable_percent);
	}

	file=fopen(fileName, "w");
	if (file==NULL)  {
		sprintf(errMsg,"Write_Lifetime_Distro(): Error opening file %s",fileName);	// Always test file open
		Error(errMsg);
	}

	fprintf(file,"%s\n",fileName);
	fprintf(file,"Frames	Time By t-alpha	Raw Count (total samples %d)	Norm Count	Total Decay	Norm Decay	Parts	Norm Parts\n",calc_nosamples2);
	for (i=0; i<useable_frames; i++) {
		fprintf(file,"%d	%.12lg	%d	%.12lg	%d	%.12lg",i+1,(i+1)*FRAMETSTEP/talpha,calc_histo[i],(double)calc_histo[i]/(calc_nosamples2),total_histo[i],norm_histo[i]/norm_histo[0]);
		if (i<useable_percent) fprintf(file,"	%d	%.12lg\n",percent[i],norm_percent[i]);
		else fprintf(file,"\n");
	}

	fclose(file);
	for (i=0; i<FRAMES; i++) free(in_2[i]);
	free(in_2);
	free(used_clus);
	free(total_histo);
	free(percent);
	free(norm_histo);
	free(norm_percent);
}

void Calc_Region_Cluster_Lifetime(int calc_nosamples, int f, int cluster_size, int *regions_parts) {
	int i, j, k, l, break_out, clust_lifetime;
	int temp_lives[10000];

	for (i=0; i<10000; i++) temp_lives[i]=0;

	l=0;

	for (i=0; i<calc_nosamples; i++) {
		if (f>=the_clusts[i][0] && f<=the_clusts[i][1]) {
			clust_lifetime=the_clusts[i][1]-the_clusts[i][0]+1;
			break_out=0;
			for (j=0; j<cluster_size; j++) {
				for (k=0; k<regions[f][no_regions[f]]; k++) {
					if (regions_parts[k]==the_clusts[i][j+2]) break_out=1;
					if (break_out==1) break;
				}
				if (break_out==1) break;
			}
			if (break_out==1) {
				if (l>10000) Error("Calc_Region_Cluster_Lifetime(): temp_lives 10000 array size not big enough");
				temp_lives[l]=clust_lifetime;
				regions_mean_cluster_lifetime[f][no_regions[f]]+=(double)clust_lifetime;
				l++;
			}

		}
	}

	regions_mean_cluster_lifetime[f][no_regions[f]]=regions_mean_cluster_lifetime[f][no_regions[f]]/l;

	fprintf(f_region_cluster_lifetimes,"%d",f);

	fprintf(f_region_cluster_lifetimes,"	%d	%d",regions[f][no_regions[f]],l);
	for (i=0; i<l; i++) fprintf(f_region_cluster_lifetimes,"	%d",temp_lives[i]);
	fprintf(f_region_cluster_lifetimes,"\n");
}

void Calc_Rg(int f, int *regions_parts) {
	int i, j;
	double dx, dy, dz, sep2;
	double sum;

	sum=0.0;
	for (i=0; i<regions[f][no_regions[f]]-1; i++) {
		for (j=i+1; j<regions[f][no_regions[f]]; j++) {
			dx=r[f][0][regions_parts[i]]+sidex*part_move[regions_parts[i]][0]-r[f][0][regions_parts[j]]-sidex*part_move[regions_parts[j]][0];
			dy=r[f][1][regions_parts[i]]+sidey*part_move[regions_parts[i]][1]-r[f][1][regions_parts[j]]-sidey*part_move[regions_parts[j]][1];
			dz=r[f][2][regions_parts[i]]+sidez*part_move[regions_parts[i]][2]-r[f][2][regions_parts[j]]-sidez*part_move[regions_parts[j]][2];
			sep2=dx*dx+dy*dy+dz*dz;
			sum+=sep2;
		}
	}

	if (regions[f][no_regions[f]]==1) printf("f%d region %d size %d Rg %lg %d\n",f,i,regions[f][no_regions[f]],regions_Rg[f][no_regions[f]],regions_parts[0]);

	sum=sum/(regions[f][no_regions[f]]*regions[f][no_regions[f]]);
	sum=sqrt(sum);
	regions_Rg[f][no_regions[f]]=sum;
}

void DFS_loop_noPBCs(int f, int in_part, int *regions_parts) {
	int j, k;

	regions_parts[regions_noPBCs[no_regions_noPBCs]]=in_part;
	regions_noPBCs[no_regions_noPBCs]++;
	used_part_noPBCs[in_part]=1;
	for (j=0; j<cnb[f][in_part]; j++) {
		if (used_bond_noPBCs[in_part][j]!=0) continue;
		if (used_part_noPBCs[bNums[f][in_part][j]]!=0) continue;
		used_bond_noPBCs[in_part][j]=1;
		for (k=0; k<cnb[f][bNums[f][in_part][j]]; k++) {
			if (bNums[f][bNums[f][in_part][j]][k]!=in_part) continue;
			used_bond_noPBCs[bNums[f][in_part][j]][k]=1;
			break;
		}
		DFS_loop_noPBCs(f,bNums[f][in_part][j],&regions_parts[0]);
	}
}

int DFS_noPBCs(int f, int *inc_parts) {	// return 0 if does not percolate, 1 if percolates
	int i, j, k;
	int *regions_parts, *inc_parts_sub;
	int through_x, through_y, through_z;
	double dx, dy, dz;
	char errMsg[1000];

	for (i=0; i<max_regions; i++) {
		regions_noPBCs[i]=0;
	}
	no_regions_noPBCs=0;
	regions_parts = malloc(N*sizeof(int));	if (regions_parts==NULL) { sprintf(errMsg,"DFS_noPBCs(): regions_parts[] malloc out of memory\n");	Error(errMsg); }
	inc_parts_sub = malloc(N*sizeof(int));	if (inc_parts_sub==NULL) { sprintf(errMsg,"DFS_noPBCs(): inc_parts_sub[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<N; i++) {
		regions_parts[i]=-1;
		inc_parts_sub[i]=0;
	}

	through_x=through_y=through_z=0;
	for (i=0; i<N; i++) {
		used_part_noPBCs[i]=-1;
		for (j=0; j<nB; j++) {
			used_bond_noPBCs[i][j]=-1;
		}
		if (inc_parts[i]==1) {
			used_part_noPBCs[i]=0;
			for (j=0; j<cnb[f][i]; j++) {
				if (inc_parts[bNums[f][i][j]]==1) {
					k=0;
					dx=r[f][0][i]-r[f][0][bNums[f][i][j]];
					dy=r[f][1][i]-r[f][1][bNums[f][i][j]];
					dz=r[f][2][i]-r[f][2][bNums[f][i][j]];
					if (PBCs==1) {
						if (dx<-halfSidex) {
							used_bond_noPBCs[i][j]=2;
							through_x=1;
							k=1;
						}
						else if (dx>halfSidex)   {
							used_bond_noPBCs[i][j]=2;
							through_x=1;
							k=1;
						}
						if (dy<-halfSidey) {
							used_bond_noPBCs[i][j]=2;
							through_y=1;
							k=1;
						}
						else if (dy>halfSidey)   {
							used_bond_noPBCs[i][j]=2;
							through_y=1;
							k=1;
						}
						if (dz<-halfSidez) {
							used_bond_noPBCs[i][j]=2;
							through_z=1;
							k=1;
						}
						else if (dz>halfSidez)   {
							used_bond_noPBCs[i][j]=2;
							through_z=1;
							k=1;
						}
					}
					if (k==0) used_bond_noPBCs[i][j]=0;
				}
			}
		}
	}
	//printf("f%d hi1\n",f);
	if (through_x==0 && through_y==0 && through_z==0) {
		free(regions_parts);
		free(inc_parts_sub);
		return 0;
	}
	//printf("f%d hi2\n",f);
	no_regions_noPBCs=0;
	for (i=0; i<N; i++) {
		if (used_part_noPBCs[i]!=0) continue;
		if (no_regions_noPBCs==max_regions) { sprintf(errMsg,"DFS_noPBCs(): no_regions_noPBCs %d == max_regions %d\n",no_regions_noPBCs,max_regions);	Error(errMsg); }
		DFS_loop_noPBCs(f,i,&regions_parts[0]);
		for (j=0; j<N; j++) inc_parts_sub[j]=0;
		for (j=0; j<regions_noPBCs[no_regions_noPBCs]; j++) {
			inc_parts_sub[regions_parts[j]]=1;
		}
		for (j=0; j<N; j++) {
			if (inc_parts_sub[j]!=1) continue;
			for (k=0; k<cnb[f][j]; k++) {
				if (used_bond_noPBCs[j][k]!=2) continue;
				if (inc_parts_sub[bNums[f][j][k]]==1) {
					free(regions_parts);
					free(inc_parts_sub);
					printf("f%d percolates x %d y %d z %d no_regions_noPBCs %d\n",f,through_x,through_y,through_z,no_regions_noPBCs);
					return 1;
				}
			}
		}
		no_regions_noPBCs++;
	}
	//printf("f%d x %d y %d z %d no_regions_noPBCs %d\n",f,through_x,through_y,through_z,no_regions_noPBCs);
	free(regions_parts);
	free(inc_parts_sub);
	return 0;
}

void DFS_loop(int f, int in_part, int *regions_parts) {
	int j, k;

	if (in_part==0) part_move[0][0]=part_move[0][1]=part_move[0][2]=0;
	regions_parts[regions[f][no_regions[f]]]=in_part;
	regions[f][no_regions[f]]++;
	used_part[in_part]=1;
	for (j=0; j<cnb[f][in_part]; j++) {
		if (used_bond[in_part][j]!=0) continue;
		if (used_part[bNums[f][in_part][j]]!=0) continue;
		used_bond[in_part][j]=1;
		for (k=0; k<cnb[f][bNums[f][in_part][j]]; k++) {
			if (bNums[f][bNums[f][in_part][j]][k]!=in_part) continue;
			used_bond[bNums[f][in_part][j]][k]=1;
			break;
		}
		part_move[bNums[f][in_part][j]][0]=part_move[in_part][0]+bond_type_x[in_part][j];
		part_move[bNums[f][in_part][j]][1]=part_move[in_part][1]+bond_type_y[in_part][j];
		part_move[bNums[f][in_part][j]][2]=part_move[in_part][2]+bond_type_z[in_part][j];
		DFS_loop(f,bNums[f][in_part][j],&regions_parts[0]);
	}
}

void DFS(int calc_no_samples, int f, int cluster_size, int *inc_parts) {
	int i, j;
	double dx, dy, dz;
	int through_x, through_y, through_z;
	int *regions_parts, *inc_parts_sub;
	int no_percolates;
	char errMsg[1000];
	no_percolates=0;

	regions_parts= malloc(N*sizeof(int));	if (regions_parts==NULL) { sprintf(errMsg,"DFS(): regions_parts[] malloc out of memory\n");	Error(errMsg); }
	inc_parts_sub = malloc(N*sizeof(int));	if (inc_parts_sub==NULL) { sprintf(errMsg,"DFS(): inc_parts_sub[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<N; i++) {
		inc_parts_sub[i]=0;
	}
	for (i=0; i<N; i++) {
		regions_parts[i]=-1;
		for (j=0; j<nB; j++) {
			bond_type_x[i][j]=bond_type_y[i][j]=bond_type_z[i][j]=0;
		}
	}
	through_x=through_y=through_z=0;
	for (i=0; i<N; i++) {
		used_part[i]=-1;
		for (j=0; j<nB; j++) {
			used_bond[i][j]=-1;
		}
		if (inc_parts[i]==1) {
			used_part[i]=0;
			for (j=0; j<cnb[f][i]; j++) {
				if (inc_parts[bNums[f][i][j]]==1) {
					used_bond[i][j]=0;
					dx=r[f][0][i]-r[f][0][bNums[f][i][j]];
					dy=r[f][1][i]-r[f][1][bNums[f][i][j]];
					dz=r[f][2][i]-r[f][2][bNums[f][i][j]];
					if (PBCs==1) {
						if (dx<-halfSidex) {
							through_x=1;
							bond_type_x[i][j]=-1;
						}
						else if (dx>halfSidex) {
							through_x=1;
							bond_type_x[i][j]=1;
						}
						if (dy<-halfSidey) {
							through_y=1;
							bond_type_y[i][j]=-1;
						}
						else if (dy>halfSidey) {
							through_y=1;
							bond_type_y[i][j]=1;
						}
						if (dz<-halfSidez) {
							through_z=1;
							bond_type_z[i][j]=-1;
						}
						else if (dz>halfSidez) {
							through_z=1;
							bond_type_z[i][j]=1;
						}
					}
				}
			}
		}
	}

	no_regions[f]=0;
	for (i=0; i<N; i++) {
		part_move[i][0]=part_move[i][1]=part_move[i][2]=0;
	}
	for (i=0; i<N; i++) {
		if (used_part[i]!=0) continue;
		if (no_regions[f]==max_regions) { sprintf(errMsg,"DFS(): no_regions[f] %d == max_regions %d\n",no_regions[f],max_regions);	Error(errMsg); }
		//printf("hi1\n");
		DFS_loop(f,i,&regions_parts[0]);
		//printf("hi2\n");
		for (j=0; j<N; j++) inc_parts_sub[j]=0;
		for (j=0; j<regions[f][no_regions[f]]; j++) {
			inc_parts_sub[regions_parts[j]]=1;
		}
		if (through_x!=0 || through_y!=0 || through_z!=0) {
			j=DFS_noPBCs(f,&inc_parts_sub[0]);
		}
		//printf("hi3\n");
		if (j==1) no_percolates++;
		else Calc_Rg(f,&regions_parts[0]);
		//printf("hi4\n");
		Calc_Region_Cluster_Lifetime(calc_no_samples,f,cluster_size,&regions_parts[0]);
		//printf("hi5\n");
		no_regions[f]++;
		for (j=0; j<N; j++) {
			part_move[j][0]=part_move[j][1]=part_move[j][2]=0;
		}
	}
	//printf("f%d regions %d percolating regions %d\n",f,no_regions[f],no_percolates);
	free(inc_parts_sub);
	free(regions_parts);
}

void Write_Clustered_Particles(char *fileName, int* **in) {
	// write clustered particles into an xmol frajectory
	FILE *file;
	char errMsg[1000];
	int e, f, i, write, remainder;

	file=fopen(fileName, "w");
	if (file==NULL)  {
		sprintf(errMsg,"Write_Clustered_Particles(): Error opening file %s",fileName);	// Always test file open
		Error(errMsg);
	}

	f=0;
	for (e=0;e<FRAMES;e++) {
		if (f==FRAMES) break;
		remainder=e%SAMPLEFREQ;
		if (remainder==0 && f<FRAMES && e>=STARTFROM) {
			write=1;
		}
		else write=0;
		if (write==1) {
			fprintf(file,"%d\nframe e %d f %d, clustered particles highlighted\n",N,e,f);
			for (i=0; i<N; i++) {
				if ((*in)[f][i]==1) {
					if (rtype[i]==1) {
						fprintf(file,"C");
					}
					else {
						fprintf(file,"D");
					}
				}
				else {
					if (rtype[i]==1) {
						fprintf(file,"A");
					}
					else {
						fprintf(file,"B");
					}
				}
				fprintf(file,"\n");
			}
			f++;
		}
	}

	fclose(file);
}

void Displacement_Distro(int overframes, int frame_start, int frame_finish, char *filenamepath, int* **in) {
	// writes displacement distribution of particles selected by in array at frame f with displacements between f and f+over_frames
	FILE *myfile;
	double dx, dy, dz, sep2, sep, thebins, *disp_distro;
	int disp_cnt;
	int f, i, k;
	int cnt;
	cnt=0;
    // This only works with cubic boxes at the moment
	disp_distro=malloc(nobins*sizeof(double));
	if (disp_distro==NULL) { Error("Displacement_Distro(): disp_distro[] malloc out of memory\n"); }
	for (i=0; i<nobins; i++) disp_distro[i]=0.0;

	for (f=0; f<FRAMES-overframes; f++) {
		if (f<frame_start) continue;
		if (f+overframes>=frame_finish) continue;
		for (i=0; i<N; i++) {
			if ((*in)[f][i]!=1) continue;
			dx=r[f][0][i]-r[f+overframes][0][i];
			dy=r[f][1][i]-r[f+overframes][1][i];
			dz=r[f][2][i]-r[f+overframes][2][i];
			if (dx<-halfSidex) { dx+=sidex; }
			else if (dx>halfSidex)   { dx-=sidex; }
			if (dy<-halfSidey) { dy+=sidey; }
			else if (dy>halfSidey)   { dy-=sidey; }
			if (dz<-halfSidez) { dz+=sidez; }
			else if (dz>halfSidez)   { dz-=sidez; }
			sep2=dx*dx+dy*dy+dz*dz;
			if (sep2 < range2) {
				sep=sqrt(sep2);
				k=(int)floor(sep/binWidth);
				disp_distro[k]+=1.0;
				cnt++;
			}
		}
		disp_cnt++;
	}

	myfile=fopen(filenamepath,"w");
	fprintf(myfile,"%s\n",filenamepath);
	if (myfile==NULL) Error("Displacement_Distro() : Error opening file\n");	// Always test file open
	printf("\nDisplacement_Distro(): printing displacement_distro(t_h=%d) to %s\n",overframes,filenamepath);
	fprintf(myfile,"length	Displacements over t=%lg (or %d frames)	Norm Displacements (samples %d)\n",overframes*FRAMETSTEP,overframes,cnt);

	for(i=0; i<nobins; i++) {
		thebins= (double)(i)*binWidth;
		fprintf(myfile,"%.12lg	%.12lg	%.12lg\n",thebins,disp_distro[i],disp_distro[i]/(cnt));
	}
	sep=0.0;
	for(i=0; i<nobins; i++) {
		sep+=disp_distro[i]*i*binWidth;
	}
	fprintf(myfile,"MSD(t=%d) is approx	%.12lg\n",overframes,sep/(cnt));

	fclose(myfile);
	free(disp_distro);
}

void Write_Continuously_In(int t_h, char *filename, int* **in) {
	// writes trajectory of particle labels for particles that are continuosly in clusters
	int i, j, k, all_frames;
	char output[1000], errMsg[1000];
	FILE *file;

	file=fopen(filename, "w");
	if (file==NULL)  {
		sprintf(errMsg,"Continuously_In(): Error opening file %s",output);	// Always test file open
		Error(errMsg);
	}

	for (i=0; i<FRAMES-t_h; i++) {
		fprintf(file,"%d\nframe %d of %d, particles continously in between frame f and f + %d\n",N,i,FRAMES-t_h,t_h);
		for (k=0; k<N; k++) {
			all_frames=(*in)[i][k];
			if (all_frames==0) {
				if (rtype[k]==1) {
					fprintf(file,"A\n");
				}
				else {
					fprintf(file,"B\n");
				}
			}
			else {
				for (j=i+1; j<=i+t_h; j++) {
					if ((*in)[j][k]==0) {
						if (rtype[k]==1) {
							fprintf(file,"A\n");
						}
						else {
							fprintf(file,"B\n");
						}
						break;
					}
					all_frames+=(*in)[j][k];
				}
				if (all_frames==t_h+1) {
					if (rtype[k]==1) {
						fprintf(file,"C\n");
					}
					else {
						fprintf(file,"D\n");
					}
				}
			}
		}
	}

	printf("written particles continuously in selected array from f to f + %d to %s",t_h,filename);

	fclose(file);
}

void Write_MSD_clustered_nonclustered(int start_frame, int end_frame, char *filename, int* **in) {
	// write MSD of particles that are continuously within clustered or non clustered domains over time which msd is measured
	int i, j, k;
	int all_frames;
	int *cnt_clust, *cnt_nonclust;
	double *sum_clust, *sum_nonclust;
	double dx, dy, dz;
	double sep2;
	char errMsg[1000];
	FILE *fout;

	cnt_clust = malloc((end_frame-start_frame)*sizeof(int));	if (cnt_clust==NULL) { sprintf(errMsg,"Write_MSD_clustered_nonclustered(): cnt_clust[] malloc out of memory\n");	Error(errMsg); }
	cnt_nonclust = malloc((end_frame-start_frame)*sizeof(int));	if (cnt_nonclust==NULL) { sprintf(errMsg,"Write_MSD_clustered_nonclustered(): cnt_nonclust[] malloc out of memory\n");	Error(errMsg); }
	sum_clust = malloc((end_frame-start_frame)*sizeof(double));	if (sum_clust==NULL) { sprintf(errMsg,"Write_MSD_clustered_nonclustered(): sum_clust[] malloc out of memory\n");	Error(errMsg); }
	sum_nonclust = malloc((end_frame-start_frame)*sizeof(double));	if (sum_nonclust==NULL) { sprintf(errMsg,"Write_MSD_clustered_nonclustered(): sum_nonclust[] malloc out of memory\n");	Error(errMsg); }
	for (i=0; i<end_frame-start_frame; i++) {
		cnt_clust[i]=0;
		cnt_nonclust[i]=0;
		sum_clust[i]=0.0;
		sum_nonclust[i]=0.0;
	}

	for (i=start_frame; i<end_frame-1; i++) {
		for (k=0; k<N; k++) {
			all_frames=(*in)[i][k];
			if (all_frames==0) {
				for (j=i+1; j<end_frame; j++) {
					all_frames+=(*in)[j][k];
					if (all_frames>0) break;
					dx=r[i][0][k]-r[j][0][k];
					dy=r[i][1][k]-r[j][1][k];
					dz=r[i][2][k]-r[j][2][k];
					if (PBCs==1) {
						if (dx<-halfSidex) { dx+=sidex; }
						else if (dx>halfSidex)   { dx-=sidex; }
						if (dy<-halfSidey) { dy+=sidey; }
						else if (dy>halfSidey)   { dy-=sidey; }
						if (dz<-halfSidez) { dz+=sidez; }
						else if (dz>halfSidez)   { dz-=sidez; }
					}
					sep2=dx*dx+dy*dy+dz*dz;
					sum_nonclust[j-i]+=sep2;
					cnt_nonclust[j-i]++;
				}
			}
			else {
				for (j=i+1; j<end_frame; j++) {
					all_frames+=(*in)[j][k];
					if (all_frames<j-i+1) break;
					dx=r[i][0][k]-r[j][0][k];
					dy=r[i][1][k]-r[j][1][k];
					dz=r[i][2][k]-r[j][2][k];
					if (PBCs==1) {
						if (dx<-halfSidex) { dx+=sidex; }
						else if (dx>halfSidex)   { dx-=sidex; }
						if (dy<-halfSidey) { dy+=sidey; }
						else if (dy>halfSidey)   { dy-=sidey; }
						if (dz<-halfSidez) { dz+=sidez; }
						else if (dz>halfSidez)   { dz-=sidez; }
					}
					sep2=dx*dx+dy*dy+dz*dz;
					sum_clust[j-i]+=sep2;
					cnt_clust[j-i]++;
				}
			}
		}
	}

	fout=fopen(filename, "w");
	if (fout==NULL)  {
		sprintf(errMsg,"MSD_slow_other() : Error opening file %s",filename);	// Always test file open
		Error(errMsg);
	}
	fprintf(fout,"%s\n",filename);
	fprintf(fout,"frame (using %d<=f<%d)	time	t by talpha	clustered	cnt	non_clusterd	cnt\n",start_frame,end_frame);

	for (i=1; i<end_frame-start_frame-1; i++) {
		fprintf(fout,"%d	%.12lg	%.12lg	%.12lg	%d	%.12lg	%d\n",i,(double)i*FRAMETSTEP,(double)i*FRAMETSTEP/talpha,sum_clust[i]/cnt_clust[i],cnt_clust[i],sum_nonclust[i]/cnt_nonclust[i],cnt_nonclust[i]);
	}

	free(cnt_clust);
	free(cnt_nonclust);
	free(sum_clust);
	free(sum_nonclust);
	fclose(fout);
}

void Do_Cluster(char *clustName, int clustSize, int centre_shell, int subClust, char *subClustStr) {
	char errMsg[1000];
	char output[1000];
	int calc_nosamples;

	Setup_ResetVars();
	mClust=0;
	printf("reading in cluster %s\n",clustName);
	Input_Cluster_Type(clustName, clustSize, subClust, subClustStr, centre_shell);
	printf("Read in dynamic TCC cluster %s - number of clusters is %d\n",clustName,mClust);
	if (mClust==0) {
		return;
	}

	if (doDynamicAnalysis==1 && doSplits==1) {
		printf("Finding splits between sights of cluster\n");
		Find_Splits(clustSize,mClust);
		printf("Found splits between sights of cluster\n");
	}

	sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.new_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs,doDynamicAnalysis,doSplits,clustName);
	calc_nosamples=Calc_Lifetimes(output,mClust,clustSize,centre_shell);
	printf("freeing the temporary cluster variables .... ");
	Free_temp(centre_shell);
	printf(" ... freed\n");

	if (doWriteLifetimeDistro==1) {
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.lives_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, doDynamicAnalysis,doSplits,clustName);
		printf("writing cluster lifetime histogram .... ");
		Write_Lifetime_Distro(output,clustSize,calc_nosamples,0,0,0,0,&clustered_parts);
		printf(" ... written\n");
	}

	if (doWriteIn==1 && doWriteLifetimeDistro==1) {
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.raw_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, doDynamicAnalysis,doSplits,clustName);
		printf("writing clustered particles .... ");
		Write_Clustered_Particles(output, &clustered_parts);
		printf(" ... written\n");
	}

	if (doDisplacementDistro==1 && doWriteLifetimeDistro==1) {
		if (ISNOTCUBIC!=0) {
        sprintf(errMsg,"Do_Cluster():ISNOTCUBIC!=0. Displacement distros can only be calculated for cubic boxes.");	// Always test file open
		Error(errMsg);

		}
        else{
            sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.disp%d.dispDistro_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, doDynamicAnalysis,doSplits,t_h,clustName);
            printf("writing displacement distribution over %d frames .... ",t_h);
            Displacement_Distro(t_h, percent_start, percent_end, output, &clustered_parts);
            printf(" ... written\n");
		}
	}

	if (doWriteCtsIn==1 && doWriteLifetimeDistro==1) {
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.cts%d.raw_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, doDynamicAnalysis,doSplits,t_h,clustName);
		printf("writing particles continuosly in clusters over %d frames following each frame.... ",t_h);
		Write_Continuously_In(t_h, output, &clustered_parts);
		printf(" ... written\n");
	}

	if (doWriteMSDClusNonClusCts==1 && doWriteLifetimeDistro==1) {
		sprintf(output,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.doDynAnal%d.doSplits%d.cts%d.msd_%s",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs, doDynamicAnalysis,doSplits,t_h,clustName);
		printf("writing MSD of particles continuosly in clusters over t.... ");
		Write_MSD_clustered_nonclustered(percent_start, percent_end, output, &clustered_parts);
		printf(" ... written\n");
	}


	printf("freeing processed cluster variables .... ");
	Free_the_clusts(calc_nosamples);
	printf(" ... freed\n\n");
}

//// START: main() routine
int main(int argc, char **argv) {
    FILE *rSizes;
    char errMsg[1000], other [1000];

	sprintf(fInputParamsName,"inputparameters.ini");
	Setup_ReadIniFile(fInputParamsName);

	if  (ISNOTCUBIC!=0){
		printf("reading box size data from %s\n",fBoxSizeName);
		rSizes=fopen(fBoxSizeName,"r");
		if(rSizes==NULL)  {
			sprintf(errMsg,"main() : Error opening file %s",fBoxSizeName);
			Error(errMsg);
		}
		fgets(other,1000,rSizes); //reads first line
		if  (ISNOTCUBIC==1) {
			Setup_ReadBox(rSizes);
			printf("sidex: %f, sidey: %f, sidez: %f\n", sidex,sidey,sidez);
		}
		if  (ISNOTCUBIC==3) {

			printf("======> Triclinic Box \n");
			// Setup_ReadBox(rSizes);
			// printf("sidex: %f, sidey: %f, sidez: %f, tilt %f\n", sidex,sidey,sidez,tilt);
		}
	}

    sprintf(fDynamicDatName,"dynamic.memsize.dat");
	Setup_ReadDynamicDat(fDynamicDatName);


	printf("Initiating variables...\n");
	sprintf(fBondsName,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,Vor,fc,PBCs);
	Setup_InitVars();
	printf("Initiated variables...reseting variables...\n");
	Setup_ResetVars();
	printf("reseting variables... starting xmol read\n");
	Setup_ReadXmol();
	printf("xmol file %s read in succesfully... starting bonds read\n",fXmolName);
	Setup_ReadBonds();
	printf("Bonds file %s read in succesfully\n",fXmolName);

	if (dyn_msp3!=-1) Do_Cluster("sp3", 3,0,0,"");
	if (dyn_msp3a!=-1) Do_Cluster("sp3a", 3,0,0,"");
	if (dyn_msp3b!=-1) Do_Cluster("sp3b", 4,0,0,"");
	if (dyn_msp3c!=-1) Do_Cluster("5A", 5,0,0,"");
	if (dyn_msp4!=-1) Do_Cluster("sp4", 4,0,0,"");
	if (dyn_msp4a!=-1) Do_Cluster("sp4a", 4,0,0,"");
	if (dyn_msp4b!=-1) Do_Cluster("sp4b", 5,0,0,"");
	if (dyn_m6A!=-1) Do_Cluster("6A", 6,0,0,"");
	if (dyn_msp5!=-1) Do_Cluster("sp5", 5,0,0,"");
	if (dyn_msp5a!=-1) Do_Cluster("sp5a", 5,0,0,"");
	if (dyn_msp5b!=-1) Do_Cluster("sp5b", 6,0,0,"");
	if (dyn_msp5c!=-1) Do_Cluster("7A", 7,0,0,"");
	if (dyn_m6Z!=-1) Do_Cluster("6Z", 6,0,1,"sp3c_i	sp3c_j");
	if (dyn_m7K!=-1) Do_Cluster("7K", 7,0,1,"sp3c_i	sp3c_j");
	if (dyn_m8A!=-1) Do_Cluster("8A", 8,0,1,"sp5b_i_1	sp5b_j_1	sp5c_i_2	sp5c_j_2	sp5b_i_3	sp5c_j_3	sp5b_i_4	sp5b_j_4	sp5c_i_5	sp5c_j_5	sp5b_i_6	sp5c_j_6");
	if (dyn_m8B!=-1) Do_Cluster("8B", 8,0,1,"sp5c");
	if (dyn_m8K!=-1) Do_Cluster("8K", 8,0,1,"sp3c_i	sp3c_j	sp3c_k");
	if (dyn_m9A!=-1) Do_Cluster("9A", 9,0,1,"sp4b_i	sp4b_j	sp4b_k");
	if (dyn_m9B!=-1) Do_Cluster("9B", 9,1,1,"sp5c_i	sp5c_j");
	if (dyn_m9K!=-1) Do_Cluster("9K", 9,1,1,"6A_i	6A_j");
	if (dyn_m10A!=-1) Do_Cluster("10A", 10,0,1,"sp4b_i	sp4b_j");
	if (dyn_m10B!=-1) Do_Cluster("10B", 10,1,1,"sp5c_i	sp5c_j	sp5c_k");
	if (dyn_m10K!=-1) Do_Cluster("10K", 10,1,1,"9K");
	if (dyn_m10W!=-1) Do_Cluster("10W", 10,0,1,"sp5b_i	sp5b_j	sp5b_k	sp5b_l	sp5b_m	sp5b_n");
	if (dyn_m11A!=-1) Do_Cluster("11A", 11,1,1,"6A_i	6A_j");
	if (dyn_m11B!=-1) Do_Cluster("11B", 11,1,1,"9B");
	if (dyn_m11C!=-1) Do_Cluster("11CD", 11,1,1,"sp5c_i	sp5c_j");
	if (dyn_m11E!=-1) Do_Cluster("11E", 11,0,1,"sp5c_i	sp5c_j	sp5c_k");
	if (dyn_m11F!=-1) Do_Cluster("11F", 11,0,1,"sp3c_i	sp3c_j	6A_i	6A_j");
	if (dyn_m11W!=-1) Do_Cluster("11W", 11,1,1,"10B");
	if (dyn_m12A!=-1) Do_Cluster("12A", 12,1,1,"11C");
	if (dyn_m12B!=-1) Do_Cluster("12BC", 12,1,1,"sp5c_i	sp5c_j	sp5c_k	sp5c_l	sp5c_m	sp5c_n");
	if (dyn_m12D!=-1) Do_Cluster("12D", 12,0,1,"sp5c_i	sp5c_j	sp5c_k	sp5c_l");
	if (dyn_m12E!=-1) Do_Cluster("12E", 12,0,1,"sp3c_i	sp3c_j	sp3c_k");
	if (dyn_m12K!=-1) Do_Cluster("12K", 12,1,1,"11A");
	if (dyn_m13A!=-1) Do_Cluster("13A", 13,1,0,"");
	if (dyn_m13B!=-1) Do_Cluster("13B", 13,1,1,"sp5c_i	sp5c_j");
	if (dyn_m13K!=-1) Do_Cluster("13K", 13,1,1,"11F	5A_i	5A_j");
	if (dyn_mFCC!=-1) Do_Cluster("FCC", 13,1,1,"sp3b_i	sp3b_j	sp3b_k	sp3b_l	sp3c_l");
	if (dyn_mHCP!=-1) Do_Cluster("HCP", 13,1,1,"sp3c_i	sp3c_j	sp3c_k");
	if (dyn_mBCC_9!=-1) Do_Cluster("BCC_9", 9,1,1,"sp4b_i_1	sp4b_j_1	6A_i_2	6A_j_2	sp4b_i_3	6A_j_3");
	if (dyn_mBCC_15!=-1) Do_Cluster("BCC_15", 15,1,1,"6A_i	6A_j	6A_k	6A_l	6A_m	6A_n");

	Setup_FreeVars();

	printf("FIN\n");

	return 0;
}