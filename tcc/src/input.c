#include "input.h"
#include "math.h"
#include "stdio.h"
#include "globals.h"
#include "iniparser.h"
#include "tools.h"

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
    RHO = iniparser_getdouble(ini, "run:number_density", -1);
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
    do11AcenXyz = iniparser_getboolean(ini, "output:11a", -1);
    do13AcenXyz = iniparser_getboolean(ini, "output:13a", -1);
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

    sidex=sidey=sidez=pow((double)N/RHO, 1.0/3.0);
    halfSidex=halfSidey=halfSidez=sidex/2.0;

    // print out values read from ini file
    printf("Xmol file name:%s Box file name:%s\n", fXmolName, fBoxSizeName);
    printf("ISNOTCUBIC %d\n",ISNOTCUBIC);
    printf("FRAMES %d N %d RHO %lg\n",FRAMES,N,RHO);
    printf("STARTFROM %d SAMPLEFREQ %d\n",STARTFROM,SAMPLEFREQ);
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

void xyz_parser(FILE *xyzfile) {

    char line[1000];
    char error_message[100];
    int i;
    int line_number, num_frames, max_particles, frame_particles;
    long file_offsets[1000];
    int valid_long = 0;

    line_number = num_frames = max_particles = 0;

    // Read in num particles
    while(feof(xyzfile) == 0) {
        frame_particles = get_long_from_string(line, &valid_long);
        if (valid_long != 1) {
            sprintf(error_message, "Unable to read XYZ file. Expected number of particles on line: %d", line_number);
            Error(error_message);
        }
        line_number += 1;
        if (frame_particles > max_particles) max_particles = frame_particles;
        file_offsets[num_frames] = ftell(xyzfile);
        for (i = 0; i < frame_particles+1; i++) {
            try_read_line_from_file(xyzfile);
            line_number += 1;
            if feof(xyzfile) {
                sprintf(error_message, "Unexpected end of file. Some particles are missing.");
                Error(error_message);
            }
        }
    }
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
            if (c=='A') particle_type[i]=1; // set particle_type array denoting cluster species
            else if (c=='B') particle_type[i]=2;
            else if (c=='C') particle_type[i]=1;
            else {
                sprintf(errMsg,"Setup_Readxyz(): unrecognized character of particle i %d from input frame %d\n",i,e);
                Error(errMsg);
            }
            if (PBCs == 1 && ISNOTCUBIC != 3) {
                wrap_particle_into_pbc(&tx, &ty, &tz);
            }
            x[i]=tx;    y[i]=ty;    z[i]=tz;    // set positions
            if (PRINTINFO==1) if (i==N-1) printf("f%d part%d %c %.5lg %.5lg %.5lg\n\n",f,i,c,x[i],y[i],z[i]);
        }
    }
}

void wrap_particle_into_pbc(double *tx, double *ty, double *tz) {
    // wrap particles back into the box
    if ((*tx) < -halfSidex) { (*tx) +=sidex; }
    else if ((*tx) > halfSidex)   { (*tx) -=sidex; }
    if ((*ty) < -halfSidey) { (*ty) +=sidey; }
    else if ((*ty) > halfSidey)   { (*ty) -=sidey; }
    if ((*tz) < -halfSidez) { (*tz) +=sidez; }
    else if ((*tz) > halfSidez)   { (*tz) -=sidez; }
}