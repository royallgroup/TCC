#ifndef TCC_STATS_H
#define TCC_STATS_H

void Stats_Init();	// initialize Stats routine
void Stats_FreeMem();	// free memory from stats variables
void Stats_SetA();	// Set arrays to true if the ith particle is a member of any clusters with this or a larger number of particles
void Stats_Analyse();	// output Cluster statistics to file
void count_gross_clusters();
void Accumulate_Stats();
void Stats_Report();
void Pop_Per_Frame(int f);
void Normalise_Populations();

#endif //TCC_STATS_H
