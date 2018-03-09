#include <stdio.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define M   5
#define N   5

int ccmp(const void *lhs, const void *rhs) {
    int c1 = *(int *) lhs;
    int c2 = *(int *) rhs;

    if (c1 < c2) return -1;
    else if (c2 < c1) return 1;
    else return 0;
}

int scmp(const void *lhs, const void *rhs) {
    return (int)(*(int **) lhs - *(int **) rhs);
}

int main(void) {
    int **tab;
    tab = (int **) malloc(M * sizeof(int *));

    for (int i = 0; i < M; i++) {
        tab[i] = (int *) malloc(N * sizeof(int));
    }

    srand((unsigned int) time(NULL));

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            tab[i][j] = rand()/10;
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", tab[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (size_t i = 0; i < M; i++) {
        qsort(tab[i], N, sizeof(int), ccmp);
    }

    qsort(tab, M, sizeof(int *), scmp);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", tab[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (size_t i = 0; i < M; i++) free(tab[i]);
    free(tab);

    return 0;
}