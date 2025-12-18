#ifndef __DEVICE_DATA_STRUCTURES_H
#define __DEVICE_DATA_STRUCTURES_H

#include "misc_def.h"

struct chemicals {
    unsigned int *chem_amounts;
    unsigned int *thresh_amounts;
    threshold_types *thresh_types;
} typedef chem_arr_t;

struct field {
    chem_id_t *chem_ids;
    unsigned int *deltas;
    unsigned int *start_bounds;
    unsigned int *end_bounds;
} typedef field_t;

struct sim_parameters {
    double max_time = -1.0;
    unsigned verbosity_bit_fields = 0;
} typedef sim_params_t;
 
struct output_stats {
    unsigned int steps_elapsed = 0;
    double time_elapsed = 0.0;
} typedef output_stats_t;

#endif