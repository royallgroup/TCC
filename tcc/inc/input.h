#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include "globals.h"

void read_ini_file();

void read_clusters_to_analyse();

void parse_box_file(int total_frames);

void get_NVT_box(FILE *read_box_file);

void get_box_file_offsets(FILE *read_box_file, int total_frames);

void get_box_size(int current_frame_number);

struct xyz_info parse_xyz_file();

void initialize_xyz_info(struct xyz_info* input_xyz_info);

void free_xyz_info(struct xyz_info* input_xyz_info);

void get_frame_coordinates_from_xyz(const struct xyz_info *input_xyz_info, int frame_number);

void get_coords_from_line(int frame_number, FILE *xyzfile, int particle);

void wrap_particle_into_pbc(double *tx, double *ty, double *tz);

#endif