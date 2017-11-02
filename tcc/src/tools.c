#include "tools.h"
// ======================================================================================================
// Here we collect all the functions that have some computational utility but do not play a mojor role in 
// the geometric contruction of the Topological Cluster Classification
// ======================================================================================================


void Dyn_Analyse_Free() {
    if (dyn_msp3!=-1) {
        free(dyn_hist_sp3);
        free(dyn_norm_hist_sp3);
        free(dyn_norm_corrected_hist_sp3);
        free(dyn_decay_norm_hist_sp3);
        free(dyn_decay_norm_corrected_hist_sp3);
        free(dyn_total_hist_sp3);
        free(dyn_norm_total_hist_sp3);
        free(dyn_decay_norm_total_hist_sp3);
        free(dyn_max_hist_sp3);
        free(dyn_norm_max_hist_sp3);
        free(dyn_norm_corrected_max_hist_sp3);
        free(dyn_decay_norm_max_hist_sp3);
        free(dyn_decay_norm_corrected_max_hist_sp3);
    }
    if (dyn_msp3a!=-1) {
        free(dyn_hist_sp3a);
        free(dyn_norm_hist_sp3a);
        free(dyn_norm_corrected_hist_sp3a);
        free(dyn_decay_norm_hist_sp3a);
        free(dyn_decay_norm_corrected_hist_sp3a);
        free(dyn_total_hist_sp3a);
        free(dyn_norm_total_hist_sp3a);
        free(dyn_decay_norm_total_hist_sp3a);
        free(dyn_max_hist_sp3a);
        free(dyn_norm_max_hist_sp3a);
        free(dyn_norm_corrected_max_hist_sp3a);
        free(dyn_decay_norm_max_hist_sp3a);
        free(dyn_decay_norm_corrected_max_hist_sp3a);
    }
    if (dyn_msp3b!=-1) {
        free(dyn_hist_sp3b);
        free(dyn_norm_hist_sp3b);
        free(dyn_norm_corrected_hist_sp3b);
        free(dyn_decay_norm_hist_sp3b);
        free(dyn_decay_norm_corrected_hist_sp3b);
        free(dyn_total_hist_sp3b);
        free(dyn_norm_total_hist_sp3b);
        free(dyn_decay_norm_total_hist_sp3b);
        free(dyn_max_hist_sp3b);
        free(dyn_norm_max_hist_sp3b);
        free(dyn_norm_corrected_max_hist_sp3b);
        free(dyn_decay_norm_max_hist_sp3b);
        free(dyn_decay_norm_corrected_max_hist_sp3b);
    }
    if (dyn_msp3c!=-1) {
        free(dyn_hist_sp3c);
        free(dyn_norm_hist_sp3c);
        free(dyn_norm_corrected_hist_sp3c);
        free(dyn_decay_norm_hist_sp3c);
        free(dyn_decay_norm_corrected_hist_sp3c);
        free(dyn_total_hist_sp3c);
        free(dyn_norm_total_hist_sp3c);
        free(dyn_decay_norm_total_hist_sp3c);
        free(dyn_max_hist_sp3c);
        free(dyn_norm_max_hist_sp3c);
        free(dyn_norm_corrected_max_hist_sp3c);
        free(dyn_decay_norm_max_hist_sp3c);
        free(dyn_decay_norm_corrected_max_hist_sp3c);
    }
    if (dyn_msp4!=-1) {
        free(dyn_hist_sp4);
        free(dyn_norm_hist_sp4);
        free(dyn_norm_corrected_hist_sp4);
        free(dyn_decay_norm_hist_sp4);
        free(dyn_decay_norm_corrected_hist_sp4);
        free(dyn_total_hist_sp4);
        free(dyn_norm_total_hist_sp4);
        free(dyn_decay_norm_total_hist_sp4);
        free(dyn_max_hist_sp4);
        free(dyn_norm_max_hist_sp4);
        free(dyn_norm_corrected_max_hist_sp4);
        free(dyn_decay_norm_max_hist_sp4);
        free(dyn_decay_norm_corrected_max_hist_sp4);
    }
    if (dyn_msp4a!=-1) {
        free(dyn_hist_sp4a);
        free(dyn_norm_hist_sp4a);
        free(dyn_norm_corrected_hist_sp4a);
        free(dyn_decay_norm_hist_sp4a);
        free(dyn_decay_norm_corrected_hist_sp4a);
        free(dyn_total_hist_sp4a);
        free(dyn_norm_total_hist_sp4a);
        free(dyn_decay_norm_total_hist_sp4a);
        free(dyn_max_hist_sp4a);
        free(dyn_norm_max_hist_sp4a);
        free(dyn_norm_corrected_max_hist_sp4a);
        free(dyn_decay_norm_max_hist_sp4a);
        free(dyn_decay_norm_corrected_max_hist_sp4a);
    }
    if (dyn_msp4b!=-1) {
        free(dyn_hist_sp4b);
        free(dyn_norm_hist_sp4b);
        free(dyn_norm_corrected_hist_sp4b);
        free(dyn_decay_norm_hist_sp4b);
        free(dyn_decay_norm_corrected_hist_sp4b);
        free(dyn_total_hist_sp4b);
        free(dyn_norm_total_hist_sp4b);
        free(dyn_decay_norm_total_hist_sp4b);
        free(dyn_max_hist_sp4b);
        free(dyn_norm_max_hist_sp4b);
        free(dyn_norm_corrected_max_hist_sp4b);
        free(dyn_decay_norm_max_hist_sp4b);
        free(dyn_decay_norm_corrected_max_hist_sp4b);
    }
    if (dyn_m6A!=-1) {
        free(dyn_hist_6A);
        free(dyn_norm_hist_6A);
        free(dyn_norm_corrected_hist_6A);
        free(dyn_decay_norm_hist_6A);
        free(dyn_decay_norm_corrected_hist_6A);
        free(dyn_total_hist_6A);
        free(dyn_norm_total_hist_6A);
        free(dyn_decay_norm_total_hist_6A);
        free(dyn_max_hist_6A);
        free(dyn_norm_max_hist_6A);
        free(dyn_norm_corrected_max_hist_6A);
        free(dyn_decay_norm_max_hist_6A);
        free(dyn_decay_norm_corrected_max_hist_6A);
    }
    if (dyn_msp5!=-1) {
        free(dyn_hist_sp5);
        free(dyn_norm_hist_sp5);
        free(dyn_norm_corrected_hist_sp5);
        free(dyn_decay_norm_hist_sp5);
        free(dyn_decay_norm_corrected_hist_sp5);
        free(dyn_total_hist_sp5);
        free(dyn_norm_total_hist_sp5);
        free(dyn_decay_norm_total_hist_sp5);
        free(dyn_max_hist_sp5);
        free(dyn_norm_max_hist_sp5);
        free(dyn_norm_corrected_max_hist_sp5);
        free(dyn_decay_norm_max_hist_sp5);
        free(dyn_decay_norm_corrected_max_hist_sp5);
    }
    if (dyn_msp5a!=-1) {
        free(dyn_hist_sp5a);
        free(dyn_norm_hist_sp5a);
        free(dyn_norm_corrected_hist_sp5a);
        free(dyn_decay_norm_hist_sp5a);
        free(dyn_decay_norm_corrected_hist_sp5a);
        free(dyn_total_hist_sp5a);
        free(dyn_norm_total_hist_sp5a);
        free(dyn_decay_norm_total_hist_sp5a);
        free(dyn_max_hist_sp5a);
        free(dyn_norm_max_hist_sp5a);
        free(dyn_norm_corrected_max_hist_sp5a);
        free(dyn_decay_norm_max_hist_sp5a);
        free(dyn_decay_norm_corrected_max_hist_sp5a);
    }
    if (dyn_msp5b!=-1) {
        free(dyn_hist_sp5b);
        free(dyn_norm_hist_sp5b);
        free(dyn_norm_corrected_hist_sp5b);
        free(dyn_decay_norm_hist_sp5b);
        free(dyn_decay_norm_corrected_hist_sp5b);
        free(dyn_total_hist_sp5b);
        free(dyn_norm_total_hist_sp5b);
        free(dyn_decay_norm_total_hist_sp5b);
        free(dyn_max_hist_sp5b);
        free(dyn_norm_max_hist_sp5b);
        free(dyn_norm_corrected_max_hist_sp5b);
        free(dyn_decay_norm_max_hist_sp5b);
        free(dyn_decay_norm_corrected_max_hist_sp5b);
    }
    if (dyn_msp5c!=-1) {
        free(dyn_hist_sp5c);
        free(dyn_norm_hist_sp5c);
        free(dyn_norm_corrected_hist_sp5c);
        free(dyn_decay_norm_hist_sp5c);
        free(dyn_decay_norm_corrected_hist_sp5c);
        free(dyn_total_hist_sp5c);
        free(dyn_norm_total_hist_sp5c);
        free(dyn_decay_norm_total_hist_sp5c);
        free(dyn_max_hist_sp5c);
        free(dyn_norm_max_hist_sp5c);
        free(dyn_norm_corrected_max_hist_sp5c);
        free(dyn_decay_norm_max_hist_sp5c);
        free(dyn_decay_norm_corrected_max_hist_sp5c);
    }
    if (dyn_m6Z!=-1) {
        free(dyn_hist_6Z);
        free(dyn_norm_hist_6Z);
        free(dyn_norm_corrected_hist_6Z);
        free(dyn_decay_norm_hist_6Z);
        free(dyn_decay_norm_corrected_hist_6Z);
        free(dyn_total_hist_6Z);
        free(dyn_norm_total_hist_6Z);
        free(dyn_decay_norm_total_hist_6Z);
        free(dyn_max_hist_6Z);
        free(dyn_norm_max_hist_6Z);
        free(dyn_norm_corrected_max_hist_6Z);
        free(dyn_decay_norm_max_hist_6Z);
        free(dyn_decay_norm_corrected_max_hist_6Z);
    }
    if (dyn_m7K!=-1) {
        free(dyn_hist_7K);
        free(dyn_norm_hist_7K);
        free(dyn_norm_corrected_hist_7K);
        free(dyn_decay_norm_hist_7K);
        free(dyn_decay_norm_corrected_hist_7K);
        free(dyn_total_hist_7K);
        free(dyn_norm_total_hist_7K);
        free(dyn_decay_norm_total_hist_7K);
        free(dyn_max_hist_7K);
        free(dyn_norm_max_hist_7K);
        free(dyn_norm_corrected_max_hist_7K);
        free(dyn_decay_norm_max_hist_7K);
        free(dyn_decay_norm_corrected_max_hist_7K);
    }
    if (dyn_m8A!=-1) {
        free(dyn_hist_8A);
        free(dyn_norm_hist_8A);
        free(dyn_norm_corrected_hist_8A);
        free(dyn_decay_norm_hist_8A);
        free(dyn_decay_norm_corrected_hist_8A);
        free(dyn_total_hist_8A);
        free(dyn_norm_total_hist_8A);
        free(dyn_decay_norm_total_hist_8A);
        free(dyn_max_hist_8A);
        free(dyn_norm_max_hist_8A);
        free(dyn_norm_corrected_max_hist_8A);
        free(dyn_decay_norm_max_hist_8A);
        free(dyn_decay_norm_corrected_max_hist_8A);
    }
    if (dyn_m8B!=-1) {
        free(dyn_hist_8B);
        free(dyn_norm_hist_8B);
        free(dyn_norm_corrected_hist_8B);
        free(dyn_decay_norm_hist_8B);
        free(dyn_decay_norm_corrected_hist_8B);
        free(dyn_total_hist_8B);
        free(dyn_norm_total_hist_8B);
        free(dyn_decay_norm_total_hist_8B);
        free(dyn_max_hist_8B);
        free(dyn_norm_max_hist_8B);
        free(dyn_norm_corrected_max_hist_8B);
        free(dyn_decay_norm_max_hist_8B);
        free(dyn_decay_norm_corrected_max_hist_8B);
    }
    if (dyn_m8K!=-1) {
        free(dyn_hist_8K);
        free(dyn_norm_hist_8K);
        free(dyn_norm_corrected_hist_8K);
        free(dyn_decay_norm_hist_8K);
        free(dyn_decay_norm_corrected_hist_8K);
        free(dyn_total_hist_8K);
        free(dyn_norm_total_hist_8K);
        free(dyn_decay_norm_total_hist_8K);
        free(dyn_max_hist_8K);
        free(dyn_norm_max_hist_8K);
        free(dyn_norm_corrected_max_hist_8K);
        free(dyn_decay_norm_max_hist_8K);
        free(dyn_decay_norm_corrected_max_hist_8K);
    }
    if (dyn_m9A!=-1) {
        free(dyn_hist_9A);
        free(dyn_norm_hist_9A);
        free(dyn_norm_corrected_hist_9A);
        free(dyn_decay_norm_hist_9A);
        free(dyn_decay_norm_corrected_hist_9A);
        free(dyn_total_hist_9A);
        free(dyn_norm_total_hist_9A);
        free(dyn_decay_norm_total_hist_9A);
        free(dyn_max_hist_9A);
        free(dyn_norm_max_hist_9A);
        free(dyn_norm_corrected_max_hist_9A);
        free(dyn_decay_norm_max_hist_9A);
        free(dyn_decay_norm_corrected_max_hist_9A);
    }
    if (dyn_m9B!=-1) {
        free(dyn_hist_9B);
        free(dyn_norm_hist_9B);
        free(dyn_norm_corrected_hist_9B);
        free(dyn_decay_norm_hist_9B);
        free(dyn_decay_norm_corrected_hist_9B);
        free(dyn_total_hist_9B);
        free(dyn_norm_total_hist_9B);
        free(dyn_decay_norm_total_hist_9B);
        free(dyn_max_hist_9B);
        free(dyn_norm_max_hist_9B);
        free(dyn_norm_corrected_max_hist_9B);
        free(dyn_decay_norm_max_hist_9B);
        free(dyn_decay_norm_corrected_max_hist_9B);
    }
    if (dyn_m9K!=-1) {
        free(dyn_hist_9K);
        free(dyn_norm_hist_9K);
        free(dyn_norm_corrected_hist_9K);
        free(dyn_decay_norm_hist_9K);
        free(dyn_decay_norm_corrected_hist_9K);
        free(dyn_total_hist_9K);
        free(dyn_norm_total_hist_9K);
        free(dyn_decay_norm_total_hist_9K);
        free(dyn_max_hist_9K);
        free(dyn_norm_max_hist_9K);
        free(dyn_norm_corrected_max_hist_9K);
        free(dyn_decay_norm_max_hist_9K);
        free(dyn_decay_norm_corrected_max_hist_9K);
    }
    if (dyn_m10A!=-1) {
        free(dyn_hist_10A);
        free(dyn_norm_hist_10A);
        free(dyn_norm_corrected_hist_10A);
        free(dyn_decay_norm_hist_10A);
        free(dyn_decay_norm_corrected_hist_10A);
        free(dyn_total_hist_10A);
        free(dyn_norm_total_hist_10A);
        free(dyn_decay_norm_total_hist_10A);
        free(dyn_max_hist_10A);
        free(dyn_norm_max_hist_10A);
        free(dyn_norm_corrected_max_hist_10A);
        free(dyn_decay_norm_max_hist_10A);
        free(dyn_decay_norm_corrected_max_hist_10A);
    }
    if (dyn_m10B!=-1) {
        free(dyn_hist_10B);
        free(dyn_norm_hist_10B);
        free(dyn_norm_corrected_hist_10B);
        free(dyn_decay_norm_hist_10B);
        free(dyn_decay_norm_corrected_hist_10B);
        free(dyn_total_hist_10B);
        free(dyn_norm_total_hist_10B);
        free(dyn_decay_norm_total_hist_10B);
        free(dyn_max_hist_10B);
        free(dyn_norm_max_hist_10B);
        free(dyn_norm_corrected_max_hist_10B);
        free(dyn_decay_norm_max_hist_10B);
        free(dyn_decay_norm_corrected_max_hist_10B);
    }
    if (dyn_m10K!=-1) {
        free(dyn_hist_10K);
        free(dyn_norm_hist_10K);
        free(dyn_norm_corrected_hist_10K);
        free(dyn_decay_norm_hist_10K);
        free(dyn_decay_norm_corrected_hist_10K);
        free(dyn_total_hist_10K);
        free(dyn_norm_total_hist_10K);
        free(dyn_decay_norm_total_hist_10K);
        free(dyn_max_hist_10K);
        free(dyn_norm_max_hist_10K);
        free(dyn_norm_corrected_max_hist_10K);
        free(dyn_decay_norm_max_hist_10K);
        free(dyn_decay_norm_corrected_max_hist_10K);
    }
    if (dyn_m10W!=-1) {
        free(dyn_hist_10W);
        free(dyn_norm_hist_10W);
        free(dyn_norm_corrected_hist_10W);
        free(dyn_decay_norm_hist_10W);
        free(dyn_decay_norm_corrected_hist_10W);
        free(dyn_total_hist_10W);
        free(dyn_norm_total_hist_10W);
        free(dyn_decay_norm_total_hist_10W);
        free(dyn_max_hist_10W);
        free(dyn_norm_max_hist_10W);
        free(dyn_norm_corrected_max_hist_10W);
        free(dyn_decay_norm_max_hist_10W);
        free(dyn_decay_norm_corrected_max_hist_10W);
    }
    if (dyn_m11A!=-1) {
        free(dyn_hist_11A);
        free(dyn_norm_hist_11A);
        free(dyn_norm_corrected_hist_11A);
        free(dyn_decay_norm_hist_11A);
        free(dyn_decay_norm_corrected_hist_11A);
        free(dyn_total_hist_11A);
        free(dyn_norm_total_hist_11A);
        free(dyn_decay_norm_total_hist_11A);
        free(dyn_max_hist_11A);
        free(dyn_norm_max_hist_11A);
        free(dyn_norm_corrected_max_hist_11A);
        free(dyn_decay_norm_max_hist_11A);
        free(dyn_decay_norm_corrected_max_hist_11A);
    }
    if (dyn_m11B!=-1) {
        free(dyn_hist_11B);
        free(dyn_norm_hist_11B);
        free(dyn_norm_corrected_hist_11B);
        free(dyn_decay_norm_hist_11B);
        free(dyn_decay_norm_corrected_hist_11B);
        free(dyn_total_hist_11B);
        free(dyn_norm_total_hist_11B);
        free(dyn_decay_norm_total_hist_11B);
        free(dyn_max_hist_11B);
        free(dyn_norm_max_hist_11B);
        free(dyn_norm_corrected_max_hist_11B);
        free(dyn_decay_norm_max_hist_11B);
        free(dyn_decay_norm_corrected_max_hist_11B);
    }
    if (dyn_m11C!=-1) {
        free(dyn_hist_11C);
        free(dyn_norm_hist_11C);
        free(dyn_norm_corrected_hist_11C);
        free(dyn_decay_norm_hist_11C);
        free(dyn_decay_norm_corrected_hist_11C);
        free(dyn_total_hist_11C);
        free(dyn_norm_total_hist_11C);
        free(dyn_decay_norm_total_hist_11C);
        free(dyn_max_hist_11C);
        free(dyn_norm_max_hist_11C);
        free(dyn_norm_corrected_max_hist_11C);
        free(dyn_decay_norm_max_hist_11C);
        free(dyn_decay_norm_corrected_max_hist_11C);
    }
    if (dyn_m11E!=-1) {
        free(dyn_hist_11E);
        free(dyn_norm_hist_11E);
        free(dyn_norm_corrected_hist_11E);
        free(dyn_decay_norm_hist_11E);
        free(dyn_decay_norm_corrected_hist_11E);
        free(dyn_total_hist_11E);
        free(dyn_norm_total_hist_11E);
        free(dyn_decay_norm_total_hist_11E);
        free(dyn_max_hist_11E);
        free(dyn_norm_max_hist_11E);
        free(dyn_norm_corrected_max_hist_11E);
        free(dyn_decay_norm_max_hist_11E);
        free(dyn_decay_norm_corrected_max_hist_11E);
    }
    if (dyn_m11F!=-1) {
        free(dyn_hist_11F);
        free(dyn_norm_hist_11F);
        free(dyn_norm_corrected_hist_11F);
        free(dyn_decay_norm_hist_11F);
        free(dyn_decay_norm_corrected_hist_11F);
        free(dyn_total_hist_11F);
        free(dyn_norm_total_hist_11F);
        free(dyn_decay_norm_total_hist_11F);
        free(dyn_max_hist_11F);
        free(dyn_norm_max_hist_11F);
        free(dyn_norm_corrected_max_hist_11F);
        free(dyn_decay_norm_max_hist_11F);
        free(dyn_decay_norm_corrected_max_hist_11F);
    }
    if (dyn_m11W!=-1) {
        free(dyn_hist_11W);
        free(dyn_norm_hist_11W);
        free(dyn_norm_corrected_hist_11W);
        free(dyn_decay_norm_hist_11W);
        free(dyn_decay_norm_corrected_hist_11W);
        free(dyn_total_hist_11W);
        free(dyn_norm_total_hist_11W);
        free(dyn_decay_norm_total_hist_11W);
        free(dyn_max_hist_11W);
        free(dyn_norm_max_hist_11W);
        free(dyn_norm_corrected_max_hist_11W);
        free(dyn_decay_norm_max_hist_11W);
        free(dyn_decay_norm_corrected_max_hist_11W);
    }
    if (dyn_m12A!=-1) {
        free(dyn_hist_12A);
        free(dyn_norm_hist_12A);
        free(dyn_norm_corrected_hist_12A);
        free(dyn_decay_norm_hist_12A);
        free(dyn_decay_norm_corrected_hist_12A);
        free(dyn_total_hist_12A);
        free(dyn_norm_total_hist_12A);
        free(dyn_decay_norm_total_hist_12A);
        free(dyn_max_hist_12A);
        free(dyn_norm_max_hist_12A);
        free(dyn_norm_corrected_max_hist_12A);
        free(dyn_decay_norm_max_hist_12A);
        free(dyn_decay_norm_corrected_max_hist_12A);
    }
    if (dyn_m12B!=-1) {
        free(dyn_hist_12B);
        free(dyn_norm_hist_12B);
        free(dyn_norm_corrected_hist_12B);
        free(dyn_decay_norm_hist_12B);
        free(dyn_decay_norm_corrected_hist_12B);
        free(dyn_total_hist_12B);
        free(dyn_norm_total_hist_12B);
        free(dyn_decay_norm_total_hist_12B);
        free(dyn_max_hist_12B);
        free(dyn_norm_max_hist_12B);
        free(dyn_norm_corrected_max_hist_12B);
        free(dyn_decay_norm_max_hist_12B);
        free(dyn_decay_norm_corrected_max_hist_12B);
    }
    if (dyn_m12D!=-1) {
        free(dyn_hist_12D);
        free(dyn_norm_hist_12D);
        free(dyn_norm_corrected_hist_12D);
        free(dyn_decay_norm_hist_12D);
        free(dyn_decay_norm_corrected_hist_12D);
        free(dyn_total_hist_12D);
        free(dyn_norm_total_hist_12D);
        free(dyn_decay_norm_total_hist_12D);
        free(dyn_max_hist_12D);
        free(dyn_norm_max_hist_12D);
        free(dyn_norm_corrected_max_hist_12D);
        free(dyn_decay_norm_max_hist_12D);
        free(dyn_decay_norm_corrected_max_hist_12D);
    }
    if (dyn_m12E!=-1) {
        free(dyn_hist_12E);
        free(dyn_norm_hist_12E);
        free(dyn_norm_corrected_hist_12E);
        free(dyn_decay_norm_hist_12E);
        free(dyn_decay_norm_corrected_hist_12E);
        free(dyn_total_hist_12E);
        free(dyn_norm_total_hist_12E);
        free(dyn_decay_norm_total_hist_12E);
        free(dyn_max_hist_12E);
        free(dyn_norm_max_hist_12E);
        free(dyn_norm_corrected_max_hist_12E);
        free(dyn_decay_norm_max_hist_12E);
        free(dyn_decay_norm_corrected_max_hist_12E);
    }
    if (dyn_m12K!=-1) {
        free(dyn_hist_12K);
        free(dyn_norm_hist_12K);
        free(dyn_norm_corrected_hist_12K);
        free(dyn_decay_norm_hist_12K);
        free(dyn_decay_norm_corrected_hist_12K);
        free(dyn_total_hist_12K);
        free(dyn_norm_total_hist_12K);
        free(dyn_decay_norm_total_hist_12K);
        free(dyn_max_hist_12K);
        free(dyn_norm_max_hist_12K);
        free(dyn_norm_corrected_max_hist_12K);
        free(dyn_decay_norm_max_hist_12K);
        free(dyn_decay_norm_corrected_max_hist_12K);
    }
    if (dyn_m13A!=-1) {
        free(dyn_hist_13A);
        free(dyn_norm_hist_13A);
        free(dyn_norm_corrected_hist_13A);
        free(dyn_decay_norm_hist_13A);
        free(dyn_decay_norm_corrected_hist_13A);
        free(dyn_total_hist_13A);
        free(dyn_norm_total_hist_13A);
        free(dyn_decay_norm_total_hist_13A);
        free(dyn_max_hist_13A);
        free(dyn_norm_max_hist_13A);
        free(dyn_norm_corrected_max_hist_13A);
        free(dyn_decay_norm_max_hist_13A);
        free(dyn_decay_norm_corrected_max_hist_13A);
    }
    if (dyn_m13B!=-1) {
        free(dyn_hist_13B);
        free(dyn_norm_hist_13B);
        free(dyn_norm_corrected_hist_13B);
        free(dyn_decay_norm_hist_13B);
        free(dyn_decay_norm_corrected_hist_13B);
        free(dyn_total_hist_13B);
        free(dyn_norm_total_hist_13B);
        free(dyn_decay_norm_total_hist_13B);
        free(dyn_max_hist_13B);
        free(dyn_norm_max_hist_13B);
        free(dyn_norm_corrected_max_hist_13B);
        free(dyn_decay_norm_max_hist_13B);
        free(dyn_decay_norm_corrected_max_hist_13B);
    }
    if (dyn_m13K!=-1) {
        free(dyn_hist_13K);
        free(dyn_norm_hist_13K);
        free(dyn_norm_corrected_hist_13K);
        free(dyn_decay_norm_hist_13K);
        free(dyn_decay_norm_corrected_hist_13K);
        free(dyn_total_hist_13K);
        free(dyn_norm_total_hist_13K);
        free(dyn_decay_norm_total_hist_13K);
        free(dyn_max_hist_13K);
        free(dyn_norm_max_hist_13K);
        free(dyn_norm_corrected_max_hist_13K);
        free(dyn_decay_norm_max_hist_13K);
        free(dyn_decay_norm_corrected_max_hist_13K);
    }
    if (dyn_mFCC!=-1) {
        free(dyn_hist_FCC);
        free(dyn_norm_hist_FCC);
        free(dyn_norm_corrected_hist_FCC);
        free(dyn_decay_norm_hist_FCC);
        free(dyn_decay_norm_corrected_hist_FCC);
        free(dyn_total_hist_FCC);
        free(dyn_norm_total_hist_FCC);
        free(dyn_decay_norm_total_hist_FCC);
        free(dyn_max_hist_FCC);
        free(dyn_norm_max_hist_FCC);
        free(dyn_norm_corrected_max_hist_FCC);
        free(dyn_decay_norm_max_hist_FCC);
        free(dyn_decay_norm_corrected_max_hist_FCC);
    }
    if (dyn_mHCP!=-1) {
        free(dyn_hist_HCP);
        free(dyn_norm_hist_HCP);
        free(dyn_norm_corrected_hist_HCP);
        free(dyn_decay_norm_hist_HCP);
        free(dyn_decay_norm_corrected_hist_HCP);
        free(dyn_total_hist_HCP);
        free(dyn_norm_total_hist_HCP);
        free(dyn_decay_norm_total_hist_HCP);
        free(dyn_max_hist_HCP);
        free(dyn_norm_max_hist_HCP);
        free(dyn_norm_corrected_max_hist_HCP);
        free(dyn_decay_norm_max_hist_HCP);
        free(dyn_decay_norm_corrected_max_hist_HCP);
    }
    if (dyn_mBCC_9!=-1) {
        free(dyn_hist_BCC_9);
        free(dyn_norm_hist_BCC_9);
        free(dyn_norm_corrected_hist_BCC_9);
        free(dyn_decay_norm_hist_BCC_9);
        free(dyn_decay_norm_corrected_hist_BCC_9);
        free(dyn_total_hist_BCC_9);
        free(dyn_norm_total_hist_BCC_9);
        free(dyn_decay_norm_total_hist_BCC_9);
        free(dyn_max_hist_BCC_9);
        free(dyn_norm_max_hist_BCC_9);
        free(dyn_norm_corrected_max_hist_BCC_9);
        free(dyn_decay_norm_max_hist_BCC_9);
        free(dyn_decay_norm_corrected_max_hist_BCC_9);
    }
    if (dyn_mBCC_15!=-1) {
        free(dyn_hist_BCC_15);
        free(dyn_norm_hist_BCC_15);
        free(dyn_norm_corrected_hist_BCC_15);
        free(dyn_decay_norm_hist_BCC_15);
        free(dyn_decay_norm_corrected_hist_BCC_15);
        free(dyn_total_hist_BCC_15);
        free(dyn_norm_total_hist_BCC_15);
        free(dyn_decay_norm_total_hist_BCC_15);
        free(dyn_max_hist_BCC_15);
        free(dyn_norm_max_hist_BCC_15);
        free(dyn_norm_corrected_max_hist_BCC_15);
        free(dyn_decay_norm_max_hist_BCC_15);
        free(dyn_decay_norm_corrected_max_hist_BCC_15);
    }
}
void Error_no_free(char *msg) { // Exit program printing error message but don't try to free any memory
    printf("\nd%d %s\n",rank,msg); 
    exit(1); 
} // Quit


void Error(char *msg) { // Exit program printing error message and trying to free any allocated memory
    printf("\nd%d %s\n",rank,msg); 
    
    Setup_FreeStaticVars();
    if (doDynamics==1) {
        Setup_FreeDynamicVars();
        Dyn_Analyse_Free();
    }
    exit(1); 
} // Quit

// resize 2D arrays
int **resize_2D_int(int **the_array, int old_row_size, int new_row_size, int new_col_size, int value) {
    int i, j;
    char errMsg[1000];
    
    the_array=realloc(the_array,new_row_size*sizeof(int *));
    if (the_array == NULL) { sprintf(errMsg,"resize_2D_int(): the_array[] out of memory old_row_size %d new_row_size %d new_col_size %d\n",old_row_size,new_row_size,new_col_size); Error_no_free(errMsg); }
    for (i=old_row_size; i<new_row_size; i++) {
        the_array[i]=malloc(new_col_size*sizeof(int));
        if (the_array == NULL) { sprintf(errMsg,"resize_2D_int(): the_array[][] out of memory\n"); Error_no_free(errMsg); }
    }
    for (i=old_row_size; i<new_row_size; i++) {
        for (j=0; j<new_col_size; j++) {
            the_array[i][j]=value;
        }
    }
    
    return the_array;
}


int *resize_1D_int(int *the_array, int old_col_size, int new_col_size) {
    int i;
    char errMsg[1000];
    
    the_array=realloc(the_array,new_col_size*sizeof(int));
    if (the_array == NULL) { sprintf(errMsg,"resize_1D_int(): the_array[] out of memory old_col_size %d new_col_size %d\n",old_col_size,new_col_size); Error_no_free(errMsg); }
    for (i=old_col_size; i<new_col_size; i++) the_array[i]=-1;
    
    return the_array;
}

double *resize_1D_double(double *the_array, int old_col_size, int new_col_size) {
    int i;
    char errMsg[1000];
    
    the_array=realloc(the_array,new_col_size*sizeof(double));
    if (the_array == NULL) { sprintf(errMsg,"resize_1D_double(): the_array[] out of memory old_col_size %d new_col_size %d\n",old_col_size,new_col_size); Error_no_free(errMsg); }
    for (i=old_col_size; i<new_col_size; i++) the_array[i]=0.0;
    
    return the_array;
}

void links() {  // sorts all the particles into cells, result given by head-of-chain and linked list arrays
    int i, ic;
    for (ic=1;ic<=ncells;ic++) head[ic]=0;
    for (i=1;i<=N;i++) {
        ic = 1 + (int)((x[i-1]+ halfSide)*invcellSide) + M*((int)((y[i-1]+halfSide)*invcellSide)) + M*M*((int)((z[i-1]+halfSide)*invcellSide));
        if (ic > ncells || ic <= 0) {
            printf("d%d i %d r_x %lg r_y %lg r_z %lg side %lg halfSide %lg ic %d ncells %d\n",rank,i-1,x[i-1],y[i-1],z[i-1],side,halfSide,ic,ncells);
            Error("links(): ic > ncells, i.e. particle coord no longer in simulation box!!\n");
        }
        llist[i]=head[ic];
        head[ic]=i;
    }
}

void Dyn_add(int *arr, int f, int clusSize, int *dyn_n, int *dyn_m, int* **dyn_l, int* **dyn_hc, int do_up, int *dyn_up, int n, int do_sub, int n_sub, int* **dyn_sub, int *sub) {
    int i, j, k;

    for (i=(*dyn_n)-1;i>-1;i--) {
        k=1;
        for (j=0;j<clusSize;j++) {
            if ((*dyn_hc)[i][j]!=arr[j]) {
                k=0;
                break;
            }
        }
        if (k==0) continue;
        if (k==1) break;
    }

    if (i==-1) {
        if ((*dyn_m)==(*dyn_n)) { 
            (*dyn_l)=resize_2D_int((*dyn_l),(*dyn_m),(*dyn_m)+incrDynamicClusters,(2*initNoLifetimes+1),-1);
            (*dyn_hc)=resize_2D_int((*dyn_hc),(*dyn_m),(*dyn_m)+incrDynamicClusters,clusSize,-1);
            if (doSubClusts==1 && do_sub==1) (*dyn_sub)=resize_2D_int((*dyn_sub),(*dyn_m),(*dyn_m)+incrDynamicClusters,n_sub,-1);
            (*dyn_m)=(*dyn_m)+incrDynamicClusters;
        }
        if (doSubClusts==1 && do_up==1) dyn_up[n]=(*dyn_n);
        for (j=0;j<clusSize;j++) (*dyn_hc)[(*dyn_n)][j]=arr[j];
        (*dyn_l)[(*dyn_n)][1]=f;
        (*dyn_l)[(*dyn_n)][2]=f;
        (*dyn_l)[(*dyn_n)][0]=1;
        if (doSubClusts==1 && do_sub==1) {
            for (j=0; j<n_sub; j++) (*dyn_sub)[(*dyn_n)][j]=sub[j];
        }
        (*dyn_n)=(*dyn_n)+1;
    }
    else {
        if (doSubClusts==1 && do_up==1) dyn_up[n]=i;
        if (doSubClusts==1 && do_sub==1) {
            for (j=0; j<n_sub; j++) {
                if ((*dyn_sub)[i][j]==-1 && sub[j]!=-1) (*dyn_sub)[i][j]=sub[j];
            }
        }
        if ((*dyn_l)[i][2*(*dyn_l)[i][0]]==f-1) (*dyn_l)[i][2*(*dyn_l)[i][0]]++;
        else {
            if ((*dyn_l)[i][0]>=initNoLifetimes) (*dyn_l)[i]=resize_1D_int((*dyn_l)[i],2*(*dyn_l)[i][0]+1,2*(*dyn_l)[i][0]+3);
            (*dyn_l)[i][2*(*dyn_l)[i][0]+1]=f;
            (*dyn_l)[i][2*(*dyn_l)[i][0]+2]=f;
            (*dyn_l)[i][0]++;
        }
    }
}

void Dyn_add_6A(int repeat, int *arr, int f, int clusSize, int *dyn_n, int *dyn_m, int* **dyn_l, int* **dyn_hc, int do_up, int *dyn_up, int n, int do_sub, int n_sub, int* **dyn_sub, int *sub) {
    int i, j, k;
    
    for (i=(*dyn_n)-1;i>-1;i--) {
        k=1;
        for (j=0;j<clusSize;j++) {
            if ((*dyn_hc)[i][j]!=arr[j]) {
                k=0;
                break;
            }
        }
        if (k==0) continue;
        if (k==1) break;
    }
    if (i==-1 && repeat==0) {
        if ((*dyn_m)==(*dyn_n)) {
            (*dyn_l)=resize_2D_int((*dyn_l),(*dyn_m),(*dyn_m)+incrDynamicClusters,(2*initNoLifetimes+1),-1);
            (*dyn_hc)=resize_2D_int((*dyn_hc),(*dyn_m),(*dyn_m)+incrDynamicClusters,clusSize,-1);
            if (doSubClusts==1 && do_sub==1) (*dyn_sub)=resize_2D_int((*dyn_sub),(*dyn_m),(*dyn_m)+incrDynamicClusters,n_sub,-1);
            (*dyn_m)=(*dyn_m)+incrDynamicClusters;
        }
        if (doSubClusts==1 && do_up==1) dyn_up[n]=(*dyn_n);
        for (j=0;j<clusSize;j++) (*dyn_hc)[(*dyn_n)][j]=arr[j];
        (*dyn_l)[(*dyn_n)][1]=f;
        (*dyn_l)[(*dyn_n)][2]=f;
        (*dyn_l)[(*dyn_n)][0]=1;
        if (doSubClusts==1 && do_sub==1) {
            for (j=0; j<n_sub; j++) (*dyn_sub)[(*dyn_n)][j]=sub[j];
        }
        (*dyn_n)=(*dyn_n)+1;
    }
    else if (i>-1 && repeat==0) {
        if (doSubClusts==1 && do_up==1) dyn_up[n]=i;
        if (doSubClusts==1 && do_sub==1) {
            for (j=0; j<n_sub; j++) {
                if ((*dyn_sub)[i][j]==-1 && sub[j]!=-1) (*dyn_sub)[i][j]=sub[j];
            }
        }
        if ((*dyn_l)[i][2*(*dyn_l)[i][0]]==f-1) (*dyn_l)[i][2*(*dyn_l)[i][0]]++;
        else {
            if ((*dyn_l)[i][0]>=initNoLifetimes) (*dyn_l)[i]=resize_1D_int((*dyn_l)[i],2*(*dyn_l)[i][0]+1,2*(*dyn_l)[i][0]+3);
            (*dyn_l)[i][2*(*dyn_l)[i][0]+1]=f;
            (*dyn_l)[i][2*(*dyn_l)[i][0]+2]=f;
            (*dyn_l)[i][0]++;
        }
    }
    else if (i>-1 && repeat==1 && doSubClusts==1 && do_up==1) dyn_up[n]=i;
}

void Dyn_add_8A(int *arr, int f, int clusSize, int *dyn_n, int *dyn_m, int* **dyn_l, int* **dyn_hc, int do_up, int *dyn_up, int n, int do_sub, int n_sub, int* **dyn_sub, int *sub) {
    int i, j, k;
    char errMsg[1000];

    for (i=(*dyn_n)-1;i>-1;i--) {
        k=1;
        for (j=0;j<clusSize;j++) {
            if ((*dyn_hc)[i][j]!=arr[j]) {
                k=0;
                break;
            }
        }
        if (k==0) continue;
        if (k==1) break;
    }

    if (sub[6]==-1 &&  sub[8]==-1 &&  sub[10]==-1) {
        if (i==-1) {
            if ((*dyn_m)==(*dyn_n)) { 
                (*dyn_l)=resize_2D_int((*dyn_l),(*dyn_m),(*dyn_m)+incrDynamicClusters,(2*initNoLifetimes+1),-1);
                (*dyn_hc)=resize_2D_int((*dyn_hc),(*dyn_m),(*dyn_m)+incrDynamicClusters,clusSize,-1);
                if (doSubClusts==1 && do_sub==1) (*dyn_sub)=resize_2D_int((*dyn_sub),(*dyn_m),(*dyn_m)+incrDynamicClusters,n_sub,-1);
                (*dyn_m)=(*dyn_m)+incrDynamicClusters;
            }
            if (doSubClusts==1 && do_up==1) dyn_up[n]=(*dyn_n);
            for (j=0;j<clusSize;j++) (*dyn_hc)[(*dyn_n)][j]=arr[j];
            (*dyn_l)[(*dyn_n)][1]=f;
            (*dyn_l)[(*dyn_n)][2]=f;
            (*dyn_l)[(*dyn_n)][0]=1;
            if (doSubClusts==1 && do_sub==1) {
                for (j=0; j<n_sub; j++) (*dyn_sub)[(*dyn_n)][j]=sub[j];
            }
            (*dyn_n)=(*dyn_n)+1;
        }
        else {
            if (doSubClusts==1 && do_up==1) dyn_up[n]=i;
            if (doSubClusts==1 && do_sub==1) {
                for (j=0; j<n_sub; j++) {
                    if ((*dyn_sub)[i][j]==-1 && sub[j]!=-1) (*dyn_sub)[i][j]=sub[j];
                }
            }
            if ((*dyn_l)[i][2*(*dyn_l)[i][0]]==f-1) (*dyn_l)[i][2*(*dyn_l)[i][0]]++;
            else {
                if ((*dyn_l)[i][0]>=initNoLifetimes) (*dyn_l)[i]=resize_1D_int((*dyn_l)[i],2*(*dyn_l)[i][0]+1,2*(*dyn_l)[i][0]+3);
                (*dyn_l)[i][2*(*dyn_l)[i][0]+1]=f;
                (*dyn_l)[i][2*(*dyn_l)[i][0]+2]=f;
                (*dyn_l)[i][0]++;
            }
        }
    }
    else {
        if (i==-1) { 
            sprintf(errMsg,"Dyn_add_8A(): Error here, re-detection of 8A in same frame doesnt appear in dyn_hc8A list");    Error(errMsg);
        }
        if (doSubClusts==1 && do_sub!=-1) {
            for (j=0; j<n_sub; j++) {
                if ((*dyn_sub)[i][j]==-1 && sub[j]!=-1) (*dyn_sub)[i][j]=sub[j];
            }
        }
    }
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BEGIN QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns YES if sort was successful, or NO if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7

int quickSort(int *arr, int elements) {
    int piv, beg[10], end[10], i, L, R;

    i=0;
    beg[0]=0; 
    end[0]=elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R) {
            piv=arr[L]; 
            if (i==10-1) return 0;
            while (L<R) {
                while (arr[R]>=piv && L<R) R--; 
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++; 
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv; 
            beg[i+1]=L+1; 
            end[i+1]=end[i]; 
            end[i++]=L; 
        }
        else i--;
    }
    return 1; 
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END QUICKSORT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


