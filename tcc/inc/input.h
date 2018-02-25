#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include "globals.h"

void Setup_ReadIniFile(char *);

void Setup_ReadBox(FILE *);

struct xyz_info parse_xyz_file(struct xyz_info input_xyz_info);

void initialize_xyz_info(struct xyz_info* input_xyz_info);

void get_xyz_frame(const struct xyz_info* input_xyz_info, int frame_number);

void get_coords_from_line(int frame_number, FILE *xyzfile, int particle);

void Setup_Readxyz(int e, int write, int f, FILE *);

void wrap_particle_into_pbc(double *tx, double *ty, double *tz);

#endif