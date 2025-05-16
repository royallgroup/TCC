#ifndef TCC_OUTPUT_H
#define TCC_OUTPUT_H
#include <stdio.h>

void Write_Raw(int f);

void Write_Raw_Particle_Types(int f, FILE *thefile, const char *sarr);

void Write_Bonds_File(int f);

void Write_Cluster_Centers_xyz(int f, int cluster_type);

void Write_Cluster_XYZ(int f);

void Write_Cluster(int f);

void Write_Cluster_Compostions(int f, int cluster_type);

void Write_Pop_Per_Frame(int f);

void write_output_files(int current_frame_number);

#endif //TCC_OUTPUT_H
