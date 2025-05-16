#ifndef RINGS_H
#define RINGS_H 

void get_basic_clusters();
void get_basic_sp4_rings(int, int, int);
void get_basic_sp5_rings(int, int, int, int);

void get_sp3_clusters(int, int, int);
void Store_sp3a(int n0, int n1, int n2);
void Store_sp3b(int n0, int n1, int n2, const int *cp);
void Store_sp3c(int n0, int n1, int n2, const int *cp);

void get_sp4_clusters(int, int, int, int);
void Store_sp4a(int n0, int n1, int n2, int n3);
void Store_sp4b(int n0, int n1, int n2, int n3, const int *cp);
void Store_sp4c(int n0, int n1, int n3, int n2, const int *cp);
void get_6A_clusters();

void get_sp5_clusters(int, int, int, int, int);
void Store_sp5a(int n0, int n1, int n2, int n3, int n4);
void Store_sp5b(int n0, int n1, int n2, int n3, int n4, const int *cp);
void Store_sp5c(int n0, int n1, int n2, int n3, int n4, const int *cp);

void add_mem_sp3b(int particle_ID);
void add_mem_sp3c(int particle_ID);
void add_mem_sp4b(int particle_ID);
void add_mem_sp4c(int particle_ID);
void add_mem_sp5b(int particle_ID);
void add_mem_sp5c(int particle_ID);


#endif
