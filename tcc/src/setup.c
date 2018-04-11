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
                fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[eleven_A_number]);
        file_pointer = open_file(output_file, "w");
        fclose(file_pointer);
    }

    if(do13AcenXyz == 1) {
        make_directory("centers_output");
        sprintf(output_file, "centers_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.%s_cen.xyz",
                fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[thirteen_A_number]);
        file_pointer = open_file(output_file, "w");
        fclose(file_pointer);
    }

    if(doWriteXYZ == 1) {
        make_directory("cluster_xyzs");
        for(cluster_number = 0; cluster_number < num_cluster_types; cluster_number++) {
            if (*do_cluster_list[cluster_number] == 1) {
                sprintf(output_file, "cluster_xyzs/%s.%s_clusts.xyz", fXmolName, cluster_names[cluster_number]);
                file_pointer = open_file(output_file, "w");
                fclose(file_pointer);
            }
        }
    }

    if(doWriteBonds == 1) {
        sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.bonds",fXmolName,rcutAA,rcutAB,rcutBB,use_voronoi_bonds,fc,PBCs);
        file_pointer=fopen(output_file, "w");
        fclose(file_pointer);
    }

    if(doWritePopPerFrame == 1) {
        sprintf(output_file,"%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.pop_per_frame",fXmolName,rcutAA,rcutAB,rcutBB,use_voronoi_bonds,fc,PBCs);
        file_pointer=fopen(output_file, "w");
        fclose(file_pointer);
    }

    if(doWriteRaw == 1) {
        make_directory("raw_output");
        for(cluster_number=0; cluster_number < num_cluster_types; cluster_number++) {
            if (*do_cluster_list[cluster_number] == 1) {
                sprintf(output_file, "raw_output/%s.rcAA%lg.rcAB%lg.rcBB%lg.Vor%d.fc%lg.PBCs%d.raw_%s",
                        fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[cluster_number]);
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
                        fXmolName, rcutAA, rcutAB, rcutBB, use_voronoi_bonds, fc, PBCs, cluster_names[cluster_number]);
                file_pointer = open_file(output_file, "w");
                fclose(file_pointer);
            }
        }
    }
}

void Initialise_Global_Variables() {
    int cluster_type;

    initNoStatic=incrStatic=1000;
    initNoClustPerPart=incrClustPerPart=1;

    mean_pop_per_frame = malloc(num_cluster_types*sizeof(double));

    x = malloc(max_particle_number*sizeof(double));   if (x==NULL) Error_no_free("Initialise_Global_Variables(): x[] malloc out of memory\n");    // positions of particles in a configuration
    y = malloc(max_particle_number*sizeof(double));   if (y==NULL) Error_no_free("Initialise_Global_Variables(): y[] malloc out of memory\n");
    z = malloc(max_particle_number*sizeof(double));   if (z==NULL) Error_no_free("Initialise_Global_Variables(): z[] malloc out of memory\n");
    particle_type=malloc(max_particle_number*sizeof(int)); if (particle_type==NULL) { Error_no_free("Initialise_Global_Variables(): particle_type[] malloc out of memory\n"); }    // type of species

    num_gross_particles = malloc(num_cluster_types*sizeof(int));
    total_clusters = malloc(num_cluster_types*sizeof(int));
    pop_per_frame = malloc(num_cluster_types*sizeof(double *));

    for(cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {

        // particle fraction of particles in each cluster in each frame
        pop_per_frame[cluster_type] = malloc(frames_to_analyse * sizeof(double));
        check_null_pointer((void *) pop_per_frame[cluster_type], "pop_per_frame");

        num_gross_particles[cluster_type] = 0;
        total_clusters[cluster_type] = 0;
    }
}

void initialise_frame_variables() {

    tiltxy = tiltxz = tiltyz = 0;

    num_bonds = malloc(max_particle_number * sizeof(int));  if (num_bonds==NULL) Error_no_free("Initialise_Global_Variables(): num_bonds[] malloc out of memory\n");    // number of "bonded" neighbours of a particle
    bNums = malloc(max_particle_number * sizeof(int *));    if (bNums==NULL) Error_no_free("Initialise_Global_Variables(): bNums[] malloc out of memory\n");    // list of bonded particles to each particle
    squared_bondlengths = malloc(max_particle_number * sizeof(double *));   if (squared_bondlengths==NULL) Error_no_free("Initialise_Global_Variables(): squared_bondlengths[] malloc out of memory\n");    // array of bond lengths


    for (int particle_number = 0; particle_number < max_particle_number; particle_number++) {
        bNums[particle_number] = malloc(max_num_bonds*sizeof(int));
        check_null_pointer((void *) bNums[particle_number], "bnums[][]");

        squared_bondlengths[particle_number] = malloc(max_num_bonds*sizeof(double));
        check_null_pointer((void *) squared_bondlengths[particle_number], "squared_bondlengths[][]");
    }

    mmem_sp3b=mmem_sp3c=mmem_sp4b=mmem_sp4c=mmem_sp5b=mmem_sp5c=initNoClustPerPart;

    // arrays for the number of clusters of each type bonded to each particle
    mem_sp3b = malloc(max_particle_number * sizeof(int *));
    if (mem_sp3b == NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp3b[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp3b[i] = malloc(mmem_sp3b * sizeof(int));
        if (mem_sp3b[i]==NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp3b[][] malloc out of memory\n");
        }
    }
    nmem_sp3b = malloc(max_particle_number * sizeof(int));
    if (nmem_sp3b == NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp3b[] malloc out of memory\n");
    }

    mem_sp3c = malloc(max_particle_number*sizeof(int *));
    if (mem_sp3c == NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp3c[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp3c[i] = malloc(mmem_sp3c * sizeof(int));
        if (mem_sp3c[i] == NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp3c[][] malloc out of memory\n");
        }
    }
    nmem_sp3c = malloc(max_particle_number * sizeof(int));
    if (nmem_sp3c == NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp3c[] malloc out of memory\n");
    }

    mem_sp4b = malloc(max_particle_number*sizeof(int *));
    if (mem_sp4b == NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp4b[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp4b[i] = malloc(mmem_sp4b*sizeof(int));
        if (mem_sp4b[i]==NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp4b[][] malloc out of memory\n");
        }
    }
    nmem_sp4b = malloc(max_particle_number*sizeof(int));
    if (nmem_sp4b==NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp4b[] malloc out of memory\n");
    }

    mem_sp4c = malloc(max_particle_number*sizeof(int *));
    if (mem_sp4c==NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp4c[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp4c[i] = malloc(mmem_sp4c * sizeof(int));
        if (mem_sp4c[i] == NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp4c[][] malloc out of memory\n");
        }
    }
    nmem_sp4c = malloc(max_particle_number * sizeof(int));
    if (nmem_sp4c == NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp4c[] malloc out of memory\n");
    }

    mem_sp5b = malloc(max_particle_number*sizeof(int *));
    if (mem_sp5b==NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp5b[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp5b[i] = malloc(mmem_sp5b*sizeof(int));
        if (mem_sp5b[i]==NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp5b[][] malloc out of memory\n");
        }
    }
    nmem_sp5b = malloc(max_particle_number*sizeof(int));
    if (nmem_sp5b==NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp5b[] malloc out of memory\n");
    }

    mem_sp5c = malloc(max_particle_number * sizeof(int *));
    if (mem_sp5c==NULL) {
        Error_no_free("Initialise_Global_Variables(): mem_sp5c[] malloc out of memory\n");
    }
    for (int i = 0; i < max_particle_number; i++) {
        mem_sp5c[i] = malloc(mmem_sp5c*sizeof(int));
        if (mem_sp5c[i] == NULL) {
            Error_no_free("Initialise_Global_Variables(): mem_sp5c[][] malloc out of memory\n");
        }
    }
    nmem_sp5c = malloc(max_particle_number*sizeof(int));
    if (nmem_sp5c==NULL) {
        Error_no_free("Initialise_Global_Variables(): nmem_sp5c[] malloc out of memory\n");
    }

    for(int cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {

        *(cluster_list_width[cluster_type]) = initNoStatic;

        // Cluster storage arrays
        *(cluster_list[cluster_type]) = malloc(*(cluster_list_width[cluster_type]) * sizeof(int *));
        check_null_pointer((void *) *(cluster_list[cluster_type]), cluster_names[cluster_type]);

        for (int i = 0; i < *(cluster_list_width[cluster_type]); ++i) {
            (*(cluster_list[cluster_type]))[i] = malloc(cluster_size[cluster_type] * sizeof(int));
            check_null_pointer((void *) (*(cluster_list[cluster_type]))[i], cluster_names[cluster_type]);
        }

        // character arrays listing what type each particle is when found in a cluster
        *(raw_list[cluster_type]) = malloc(max_particle_number * sizeof(char));
        check_null_pointer((void *) *(raw_list[cluster_type]), cluster_names[cluster_type]);
    }
}

void check_null_pointer(void *pointer, char *pointer_name) {
    char errMsg[200];
    if (pointer == NULL) {
            sprintf(errMsg,"Initialise_Global_Variables(): %s malloc out of memory\n", pointer_name);
            Error_no_free(errMsg);
        }
}

void Reset_Frame_Variables() { // Reset static variables in each frame
    int cluster_type;

    for(cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        *num_cluster_list[cluster_type] = 0;
    }

    memset(num_bonds, 0, particles_in_current_frame* sizeof(int));

    memset(nmem_sp3b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp3c, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp4b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp4c, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp5b, 0, sizeof(int)*max_particle_number);
    memset(nmem_sp5c, 0, sizeof(int)*max_particle_number);

    for(cluster_type=0; cluster_type < num_cluster_types; cluster_type++) {
        memset(*raw_list[cluster_type], 'C', max_particle_number*sizeof(char));
    }
}

void Free_All_Variables()  {  // Free bond detection variables
    int i;

    free(mean_pop_per_frame);

    free(particle_type);
    free(fXmolName);
    free(fBoxSizeName);
    free(x); free(y); free(z);

    for (i = 0; i < max_particle_number; i++) {
        free(bNums[i]); 
        free(squared_bondlengths[i]);
    }
    free(bNums); free(squared_bondlengths); free(num_bonds);

    for (i = 0; i < max_particle_number; i++) free(mem_sp3b[i]);
    for (i = 0; i < max_particle_number; i++) free(mem_sp3c[i]);
    for (i = 0; i < max_particle_number; i++) free(mem_sp4b[i]);
    for (i = 0; i < max_particle_number; i++) free(mem_sp4c[i]);
    for (i = 0; i < max_particle_number; i++) free(mem_sp5b[i]);
    for (i = 0; i < max_particle_number; i++) free(mem_sp5c[i]);
    
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


    for(int cluster_type = 0; cluster_type < num_cluster_types; cluster_type++) {
        for (i = 0; i < *(cluster_list_width[cluster_type]); i++) {
            free((*cluster_list[cluster_type])[i]);
        }
        free(pop_per_frame[cluster_type]);
        free(*raw_list[cluster_type]);
        free(*cluster_list[cluster_type]);

    }

    free(pop_per_frame);
    free(num_gross_particles);
    free(total_clusters);
}

void analyse_cluster_dependencies() {

    if(doBCC15 == 1) dosp4c = 1;
    if(doBCC9 == 1) dosp4b = dosp4c = 1;
    if(doFCC == 1) dosp3b = dosp3c = 1;
    if(doHCP == 1) dosp3c = 1;
    if(do13K == 1) dosp3c = dosp4c = do11F = 1;
    if(do13B == 1) dosp5c = 1;
    if(do13A == 1) do12B = 1;
    if(do12K == 1) do11A = 1;
    if(do12E == 1) dosp3c = do11F = 1;
    if(do12D == 1) dosp5c = do11E = 1;
    if(do12B == 1) dosp5c = 1;
    if(do12A == 1) do11C = 1;
    if(do11W == 1) do10B = 1;
    if(do11F == 1) dosp3c = dosp4c = 1;
    if(do11E == 1) dosp5c = do9B = 1;
    if(do11C == 1) dosp5c = 1;
    if(do11B == 1) do9B = 1;
    if(do11A == 1) dosp4c = 1;
    if(do10W == 1) dosp5b = 1;
    if(do10K == 1) do9K = 1;
    if(do10B == 1) dosp5c = do9K = 1;
    if(do10A == 1) dosp4b = 1;
    if(do9K == 1) dosp4c = 1;
    if(do9B == 1) dosp5c = 1;
    if(do9A ==1) dosp4b = 1;
    if(do8K == 1) dosp3c =1;
    if(do8B == 1) dosp5c = 1;
    if(do8A == 1) dosp5b = dosp5c = 1;
    if(do7T_a == 1 || do7T_s == 1) do6Z = 1;
    if(do7K == 1) dosp3c = 1;
    if(do6Z == 1) dosp3c = 1;
    if(do6A == 1) dosp4c = 1;
    if(dosp5c == 1) dosp5 = 1;
    if(dosp5b == 1) dosp5 = 1;
    if(dosp5a == 1) dosp5 = 1;
    if(dosp4c == 1) dosp4 = 1;
    if(dosp4b == 1) dosp4 = 1;
    if(dosp4a == 1) dosp4 = 1;
    if(dosp3c == 1) dosp3 = 1;
    if(dosp3b == 1) dosp3 = 1;
    if(dosp3a == 1) dosp3 = 1;
    if(dosp5 == 1) dosp4 = 1;
    if(dosp4 == 1) dosp3 = 1;

}

void setup_cluster_lists() {
    num_cluster_types = 0;

    while(strcmp(cluster_names[num_cluster_types], "-1") != 0) {
        num_cluster_types += 1;
    }

    if(cluster_size[num_cluster_types] != -1)
        Error_no_free("Incorrect number of elements in cluster_size array");

    if(do_cluster_list[num_cluster_types] != NULL)
        Error_no_free("Incorrect number of elements in do_cluster_list");

    if(num_cluster_list[num_cluster_types] != NULL)
        Error_no_free("Incorrect number of elements in num_cluster_list");

    if(raw_list[num_cluster_types] != NULL)
        Error_no_free("Incorrect number of elements in raw_list");

    if(cluster_list[num_cluster_types] != NULL)
        Error_no_free("Incorrect number of elements in cluster_list");

    if(cluster_list_width[num_cluster_types] != NULL)
        Error_no_free("Incorrect number of elements in cluster_list_width");
}