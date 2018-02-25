#include "input.h"
#include "math.h"
#include "stdio.h"
#include "globals.h"
#include "iniparser.h"
#include "tools.h"

void Setup_ReadIniFile(char *filename) {

    char errMsg[1000];
    dictionary  *   ini ;

    fXmolName=malloc(500*sizeof(char)); if (fXmolName==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): fXmolName[] malloc out of memory\n");   Error_no_free(errMsg); }
    fBoxSizeName=malloc(500*sizeof(char)); if (fBoxSizeName==NULL) { sprintf(errMsg,"Initialise_Global_Variables(): fBoxSizeName[] malloc out of memory\n");   Error_no_free(errMsg); }

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
    RHO = iniparser_getdouble(ini, "run:number_density", -1);
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

    // print out values read from ini file
    printf("Xmol file name:%s Box file name:%s\n", fXmolName, fBoxSizeName);
    printf("ISNOTCUBIC %d\n",ISNOTCUBIC);
    printf("FRAMES %d SAMPLEFREQ %d\n",FRAMES, SAMPLEFREQ);
    printf("rcutAA %lg rcutAB %lg rcutBB %lg\n",rcutAA,rcutAB,rcutBB);
    printf("rcutAA2 %lg rcutAB2 %lg rcutBB2 %lg\n",rcutAA2,rcutAB2,rcutBB2);
    printf("Vor %d PBCs %d fc %lg nB %d USELIST %d\n",Vor,PBCs,fc,nB,USELIST);
    printf("write bonds file %d doWriteClus %d doWriteRaw %d doWritePopPerFrame %d\n",doWriteBonds,doWriteClus,doWriteRaw,doWritePopPerFrame);
    printf("PRINTINFO %d\n", PRINTINFO);

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

struct xyz_info parse_xyz_file(struct xyz_info input_xyz_info) {

    char line[1000];
    char error_message[100];
    int i;
    int line_number;
    int valid_long = 0;
    FILE *xyzfile;

    line_number = 0;

    initialize_xyz_info(&input_xyz_info);

    xyzfile=fopen(fXmolName,"rb");    // open xmol trajecotry
    if (xyzfile==NULL)  {
        sprintf(error_message,"Error opening XYZ file %s",fXmolName);    // Always test file open
        Error_no_free(error_message);
    }

    while(feof(xyzfile) == 0) {
        // Read in num particles
        if(input_xyz_info.total_frames > 1000) {
            sprintf(error_message, "XYZ file has over 1000 frames, cannot process XYZ file");
            Error_no_free(error_message);
        }
        line[0] = '\n';
        fgets(line, 1000, xyzfile);
        if (line[0] != '\n') {
            input_xyz_info.num_particles[input_xyz_info.total_frames] = get_long_from_string(line, &valid_long);
            if (valid_long != 1) {
                sprintf(error_message, "Unable to read XYZ file. Expected number of particles on line: %d",
                        line_number);
                Error_no_free(error_message);
            }
            line_number += 1;
            input_xyz_info.frame_offsets[input_xyz_info.total_frames] = ftell(xyzfile);
            for (i = 0; i < input_xyz_info.num_particles[input_xyz_info.total_frames]+1; i++) {
                if feof(xyzfile) {
                    sprintf(error_message, "Unexpected end of file at line %d. Some particles are missing.",
                            line_number);
                    Error_no_free(error_message);
                }
                try_read_line_from_file(xyzfile);
                line_number += 1;

            }
            input_xyz_info.total_frames += 1;
        }
    }
    fclose(xyzfile);
    return input_xyz_info;
}

void initialize_xyz_info(struct xyz_info* input_xyz_info) {
    (*input_xyz_info).total_frames = 0;
    (*input_xyz_info).data_width = 1000;
    (*input_xyz_info).num_particles = malloc((*input_xyz_info).data_width * sizeof(long));
    (*input_xyz_info).frame_offsets = malloc((*input_xyz_info).data_width * sizeof(long));
}

void get_xyz_frame(const struct xyz_info* input_xyz_info, int frame_number) {

    FILE *xyzfile;
    char error_message[1000];
    char line[1000];
    int particle;

    xyzfile=fopen(fXmolName,"rb");
    if (xyzfile==NULL)  {
        sprintf(error_message,"Error opening XYZ file %s",fXmolName);
        Error(error_message);
    }

    fseek(xyzfile, input_xyz_info->frame_offsets[frame_number], SEEK_SET);
    fgets(line, 1000, xyzfile); // comment line
    for(particle=0; particle<input_xyz_info->num_particles[frame_number]; particle++) {
        get_coords_from_line(frame_number, xyzfile, particle);
    }

}

void get_coords_from_line(int frame_number, FILE *xyzfile, int particle) {
    char line[1000], error_message[1000];
    int valid_double = 0;
    double temp_coord[3];
    char * word;
    int dimension;

    fgets(line, 1000, xyzfile);
    word = strtok (line," \t");
    if (strcmp(word, "A") == 0 || strcmp(word, "a") == 0 || strcmp(word, "1") == 0) {
        particle_type[particle] = 1;
    }
    else{
        particle_type[particle] = 2;
    }

    for(dimension=0; dimension<3; dimension++) {
        word = strtok(NULL, " \t");
        temp_coord[dimension] = get_double_from_string(word, &valid_double);
        if (valid_double != 1) {
            sprintf(error_message, "Unable to read XYZ file. Expected particle coordinate in frame number"
                    " %d particle number %d. Actually got %s", frame_number, particle, word);
            Error(error_message);
        }
    }

    if (PBCs == 1 && ISNOTCUBIC != 3) {
        wrap_particle_into_pbc(&temp_coord[0], &temp_coord[1], &temp_coord[2]);
    }
    x[particle] = temp_coord[0];
    y[particle] = temp_coord[1];
    z[particle] = temp_coord[2];
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