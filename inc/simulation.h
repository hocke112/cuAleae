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
 * This file contains verbosity bit definitions, simulation error codes, and the function prototype for the simulation master. Verbosity
 * bit fields originated from aleae.h of Aleae's source code.
*/

#ifndef __SIMULATION_HPP
#define __SIMULATION_HPP

#include "device_data_structures.h"

// bit-field for verbosity
#define PRINT_TRIALS     1
#define PRINT_TERMINAL   2
#define PRINT_TRACE      4
#define PRINT_STATES     8

enum simulation_err_t {
    SIMULATION_SUCCESS   = 0,
    INVALID_REACT_CHOSEN = 1,
    REAC_BOUND_ARR_ERR   = 2,
    PROD_BOUND_ARR_ERR   = 3,
    THRESH_CODE_ERR      = 4
};

extern "C" simulation_err_t simulation_master(unsigned int *post_trial_chem_amounts, unsigned int *total_count_triggered_threshs_h,
                                                crn_t crn_h, sim_params_t sim_params, output_stats_t *out_stats,
                                                unsigned int trial_num);

#endif