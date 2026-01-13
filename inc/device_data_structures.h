/*
 * Copyright (c) 2025-2026 Jeremiah Hockett, Chiemeka Nwakama, Aditya Parida, and Cheo Cedillo
 *
 * This file is part of cuAleae.
 *
 * cuAleae is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * cuAleae is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with cuAleae. If not, see <https://www.gnu.org/licenses/>.
 *
 *
 * This file contains definitions of structs used to hold device-side arrays.
*/

#ifndef __DEVICE_DATA_STRUCTURES_H
#define __DEVICE_DATA_STRUCTURES_H

#include "misc_def.h"

struct chemical_arrays {
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

struct crn_cu {
    chem_arr_t chem_arrays;
    field_t reactants;
    field_t products;
    float *rates;
    size_t num_chems = 0;
    size_t num_reactions = 0;
    size_t num_reactants = 0;
    size_t num_products = 0;
} typedef crn_t;

struct sim_parameters {
    double max_time = -1.0;
    int max_steps = -1;
    unsigned verbosity_bit_fields = 0;
} typedef sim_params_t;
 
struct output_stats {
    unsigned int steps_elapsed = 0;
    double time_elapsed = 0.0;
} typedef output_stats_t;

#endif