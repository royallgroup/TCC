#ifndef SETUP_H
#define SETUP_H

void Setup_Output_Files();

void initialise_run_variables();

void initialise_frame_variables(int num_particles);

void check_null_pointer(void *pointer, char *pointer_name);

void Reset_Frame_Variables();

void free_run_variables();

void free_frame_variables(int num_particles);

void analyse_cluster_dependencies();

void setup_cluster_lists();

#endif