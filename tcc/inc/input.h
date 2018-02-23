#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

void Setup_ReadIniFile(char *);
void Setup_ReadBox(FILE *);
void xyz_parser(FILE *xyzfile);
void Setup_Readxyz(int e, int write, int f, FILE *);  // output single float valued arrays in gopenmol xyz format
void wrap_particle_into_pbc(double *tx, double *ty, double *tz);

#endif