#ifndef __SIMULATION_HPP
#define __SIMULATION_HPP

#include "device_data_structures.h"

// bit-field for verbosity
#define PRINT_TRIALS     1
#define PRINT_TERMINAL   2
#define PRINT_TRACE      4
#define PRINT_STATES     8

extern "C" void simulation_master(unsigned int *post_trial_chem_amounts,
                                    chem_arr_t chem_arrays_h, field_t reactants_h, field_t products_h, float *rates_h,
                                    unsigned int num_chems, unsigned int num_reactions,
                                    unsigned int num_reactants, unsigned int num_products,
                                    sim_params_t sim_params, output_stats_t *out_stats,
                                    unsigned int max_steps, unsigned int *triggered_thresh, bool *within_threshold,
                                    unsigned int trial_num, int *err);

#endif