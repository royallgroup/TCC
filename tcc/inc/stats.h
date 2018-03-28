#ifndef TCC_STATS_H
#define TCC_STATS_H

#include "stdio.h"

void Stats_Init();

void Stats_FreeMem();

void count_number_of_clusters();

void Stats_Report();

FILE *open_stats_report_file(const char *output_name);

void stats_report_bonds(FILE *output_file);

void stats_report_clusters(FILE *output_file);

void stats_report_title(const char *output_name, FILE *output_file);

void Pop_Per_Frame(int f);

void count_mean_pop_per_frame(int frames_analysed);

#endif
