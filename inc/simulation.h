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

// extern "C" void simulation_master(unsigned int *post_trial_chem_amounts,
//                                     chem_arr_t chem_arrays_h, field_t reactants_h, field_t products_h, float *rates_h,
//                                     unsigned int num_chems, unsigned int num_reactions,
//                                     unsigned int num_reactants, unsigned int num_products,
//                                     sim_params_t sim_params, output_stats_t *out_stats,
//                                     unsigned int *triggered_thresh, bool *within_threshold,
//                                     unsigned int trial_num, int *err);

extern "C" void simulation_master(unsigned int *post_trial_chem_amounts, unsigned int *total_triggered_threshs_h,
                                        simulation_err_t *err, crn_t crn_h, sim_params_t sim_params, output_stats_t *out_stats,
                                        bool *within_threshold, unsigned int trial_num);

#endif