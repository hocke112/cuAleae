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
 * This file contains the CUDA C part of cuAleae. Implementation was based on aleae_stoch.cc of Aleae's source code.
 *
 * Note: Jeremiah's implementation of MP4-1 and MP4-2 for EE 5351 was used for this file.
*/

#include "device_data_structures.h"
#include "misc_def.h"
#include "simulation.h"
#include <stdio.h>
#include <cfloat>

#define BLOCK_SIZE 16
#define TILE_SIZE 1024

/*
 * A device function to calculate n choose r, which is taken from a Geeks-for-Geeks article on
 * implementations of nCr (https://www.geeksforgeeks.org/dsa/program-calculate-value-ncr/).
 *
 * @param n: an integer for the total number of choices
 * @param r: an integer for the amount of choices
 * @return the result of n choose r
 *
*/
__device__ float nCr(unsigned int n, unsigned int r){
    float sum = 1.0f;

    // Calculate the n choose r using the binomial coefficient formula
    for (unsigned int i = 1; i <= r; i++){
        sum = sum * (n - r + i) / i;
    }

    return sum;
}

/*
 * A kernel function to find the propensity of each reaction, which is based on equation 1 from the paper
 * Maximum Probability Reaction Sequences in Stochastic Chemical Kinetic Systems by Maryam Salehi and
 * Theodore J Perkins.
 *
 * @param propensities: an output array of all the propensities
 * @param chem_arrays: a struct containing arrays for the chemical amounts, threshold amounts, and threshold types
 * @param reactants: a struct containing the chemicals and coefficents found in the reactant field of a chemical equation,
 *                   as well as bound arrays to access them
 * @param rates: an input array containing the rate of each chemical reaction
 * @param num_reactions: the count of all reactions in the CRN
*/
__global__ void calculate_propensity_kernel(float* propensities, chem_arr_t chem_arrays, field_t reactants,
                                            float *rates, unsigned int num_reactions) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < num_reactions) {
        float local_propensity = rates[i];
        unsigned int start_bound = reactants.start_bounds[i];
        unsigned int end_bound = reactants.end_bounds[i];

        for (int it = start_bound; it < end_bound; ++it) {
            unsigned int x = chem_arrays.chem_amounts[reactants.chem_ids[it]];
            unsigned int c = reactants.deltas[it];

            if (x < c) {
                local_propensity = 0.0f;
                break;
            } else {
                local_propensity *= nCr(x, c);
            }
        }

        propensities[i] = max(local_propensity, 0.0f);
    }
    __syncthreads();
}

/*
 * A kernel function to find the partial sum of an arbitrarily-sized array of floats, based on Jeremiah Hockett's
 * integer prefix sum implementation of MP4-2 for EE 5351.
 *
 * @param out_array: an output array of partial sums
 * @param in_array: an input array of floats to find the prefix sum
 * @param block_sums: the prefix sum of each block
 * @param num_elements: the count of all elements in the output and input arrays
*/
__global__ void prescanArrayKernel(float *out_array, float *in_array, float *block_sums, int num_elements)
{
	__shared__ float scan_array[2*BLOCK_SIZE];

	unsigned int t = threadIdx.x;											// Calculate initial variables
	unsigned int start = 2 * blockDim.x * blockIdx.x;

	unsigned int in_data_index0 = start + t;
	unsigned int in_data_index1 = start + blockDim.x + t;

	if (in_data_index0 < num_elements)										// Load in input data to shared memory
		scan_array[t] = in_array[in_data_index0];
	else
		scan_array[t] = 0.0f;

	if (in_data_index1 < num_elements)
		scan_array[blockDim.x + t] = in_array[in_data_index1];
	else
		scan_array[blockDim.x + t] = 0.0f;

	__syncthreads();

	unsigned int stride = 1;
	while (stride <= BLOCK_SIZE) {
		unsigned int index = (t+1)*stride*2 - 1;
		if (index < 2*BLOCK_SIZE)
			scan_array[index] += scan_array[index-stride];

		stride = stride*2;
		__syncthreads();
	}

    if (t == 0) {															// Store block sums
		block_sums[blockIdx.x] = scan_array[2 * blockDim.x - 1];
	}

	stride = BLOCK_SIZE/2;
	while (stride > 0) {
		unsigned int index = (t+1)*stride*2 - 1;
		if((index+stride) < 2*BLOCK_SIZE)
			scan_array[index+stride] += scan_array[index];

		stride = stride/2;
		__syncthreads();
	}

	if (in_data_index0 < num_elements)										// Store data from shared memory to output
		out_array[in_data_index0] = scan_array[t];

	if (in_data_index1 < num_elements)
		out_array[in_data_index1] = scan_array[blockDim.x + t];
}

/*
 * A kernel function to find the partial sum of a block sums, based on Jeremiah Hockett's
 * integer prefix sum implementation of MP4-2 for EE 5351.
 *
 * @param scanned_block_sums: an output array of the scanned partial sums of each block
 * @param block_sums: an array containing prefix sum of each block
 * @param num_elements: the count of all elements in block_sums
*/
__global__ void scanBlockSumsKernel(float *scanned_block_sums, float *block_sums, int num_elements)
{
	__shared__ float scan_array[2*BLOCK_SIZE];

	unsigned int t = threadIdx.x;											// Calculate initial variables
	unsigned int start = 2 * blockDim.x * blockIdx.x;

	unsigned int in_data_index0 = start + t;
	unsigned int in_data_index1 = start + blockDim.x + t;

	if (in_data_index0 < num_elements)										// Load in unscanned block sums to shared memory
		scan_array[t] = block_sums[in_data_index0];
	else
		scan_array[t] = 0.0f;

	if (in_data_index1 < num_elements)
		scan_array[blockDim.x + t] = block_sums[in_data_index1];
	else
		scan_array[blockDim.x + t] = 0.0f;

	__syncthreads();

	unsigned int stride = 1;
	while (stride <= BLOCK_SIZE) {
		unsigned int index = (t+1)*stride*2 - 1;
		if (index < 2*BLOCK_SIZE)
			scan_array[index] += scan_array[index-stride];

		stride = stride*2;
		__syncthreads();
	}

    if (t == 0) {															// Store scannned block sums
		scanned_block_sums[blockIdx.x] = scan_array[2 * blockDim.x - 1];
	}

	stride = BLOCK_SIZE/2;
	while (stride > 0) {
		unsigned int index = (t+1)*stride*2 - 1;
		if((index+stride) < 2*BLOCK_SIZE)
			scan_array[index+stride] += scan_array[index];

		stride = stride/2;
		__syncthreads();
	}

	if (in_data_index0 < num_elements)										// Store data from shared memory to output
		block_sums[in_data_index0] = scan_array[t];

	if (in_data_index1 < num_elements)
		block_sums[in_data_index1] = scan_array[blockDim.x + t];
}

/*
 * A kernel function to add the block sum of block i to the scanned elements of block i + 1, based on Jeremiah Hockett's
 * integer prefix sum implementation of MP4-2 for EE 5351.
 *
 * @param out_array: an output array of partial sums
 * @param block_sums: an input array containing prefix sum of each block
 * @param num_out_elems: the count of all elements in out_array
 * @param num_block_sum_elems: the count of all elements in block_sums
*/
__global__ void addScanBlocksKernel(float *out_array, float *block_sums, unsigned int num_out_elems, unsigned int num_block_sum_elems) {
	// block sizes for this kernel function has to be 2 * BLOCK_SIZE in order for it to work.
	__shared__ float block_sum;

	if (threadIdx.x == 0 && blockIdx.x < num_block_sum_elems) {				// Stash block sum into shared memory
		block_sum = block_sums[blockIdx.x];
	}

	__syncthreads();
	unsigned int out_index = blockDim.x * (blockIdx.x + 1) + threadIdx.x;	// Find and store total of current output_value and block sum
	if (out_index < num_out_elems) {
		out_array[out_index] += block_sum;
	}
}

/*
 * A kernel function to choose a reaction based on a randomly-generated number and the propensity inclusive sums. This is needed
 * for choose_reaction() to work.
 *
 * @param candidate_reactions: an array of unsigned integers representing a reaction
 * @param propensity_inclusive_sums: an input array containing inclusive prefix sum of the propensities
 * @param propensity_sum: the sum of all propensities
 * @param r: a randomly-generated number
 * @param num_reactions: the count of all reactions in the CRN
*/
__global__ void find_reaction_candidates(unsigned int *candidate_reactions, float *propensity_inclusive_sums, float propensity_sum, float r, unsigned int num_reactions) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < num_reactions) {
        float probability = (float) propensity_inclusive_sums[i] / (float) propensity_sum;

        if (probability > r) {
            candidate_reactions[i] = i;
        } else {
            candidate_reactions[i] = UINT_MAX;
        }
    }
}

/*
 * A kernel function to choose a reaction via reduction.
 *
 * @param candidate_reactions: an array of unsigned integers representing a reaction
 * @param num_reactions: the count of all reactions in the CRN
*/
__global__ void choose_reaction_kernel(unsigned int *candidate_reactions, unsigned int num_reactions) {
    __shared__ unsigned int partial_min[2 * BLOCK_SIZE]; 					// shmem for each block

	unsigned int t = threadIdx.x; 											// Setup vars for indexing to global and shmem
	unsigned int start = 2 * blockDim.x * blockIdx.x;

	unsigned int in_data_index0 = start + t;
	unsigned int in_data_index1 = start + blockDim.x + t;

	if (in_data_index0 < num_reactions)													// To prevent out-of-bounds errors
		partial_min[t] = candidate_reactions[in_data_index0]; 							// These errors occur when n is not a power of 2
	else
		partial_min[t] = UINT_MAX;

	if (in_data_index1 < num_reactions)
		partial_min[blockDim.x + t] = candidate_reactions[in_data_index1];
	else
		partial_min[blockDim.x + t] = UINT_MAX;

	for (unsigned int stride = blockDim.x; stride >= 1; stride >>= 1) { 	// Choose the reaction for each stride
		__syncthreads();
		if (t < stride)
			partial_min[t] = min(partial_min[t], partial_min[t + stride]);
	}

    if (threadIdx.x == 0)
	    candidate_reactions[blockIdx.x] = partial_min[0];									// Write back to global memory
}

/*
 * A kernel function to remove spent chemicals. This is mainly created to avoid a cudaMemcpy call of the
 * chemical amount array for each iteration of the while loop in simulation_master.
 *
 * @param chem_arrays: a struct containing arrays for the chemical amounts, threshold amounts, and threshold types
 * @param reactants: a struct containing the chemicals and coefficents found in the reactant field of a chemical equation,
 *                   as well as bound arrays to access them
 * @param start_bound: the index of the first chemical of the reactants in a reaction
 * @param end_bound: the index after the last chemical of the reactants in a reaction
*/
__global__ void update_chem_amounts_using_reactants(chem_arr_t chem_arrays, field_t reactants, int start_bound, int end_bound) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x + start_bound;

    if (i < end_bound) {
        // x normally should be unique to each thread, but in case it isn't guard against race conditions with atomics
        chem_id_t x = reactants.chem_ids[i];
        atomicSub(&chem_arrays.chem_amounts[x], reactants.deltas[i]);
    }
    __syncthreads();
}

/*
 * A kernel function to remove spent chemicals. This is mainly created to avoid a cudaMemcpy call of the
 * chemical amount array for each iteration of the while loop in simulation_master.
 *
 * @param chem_arrays: a struct containing arrays for the chemical amounts, threshold amounts, and threshold types
 * @param products: a struct containing the chemicals and coefficents found in the reactant field of a chemical equation,
 *                   as well as bound arrays to access them
 * @param start_bound: the index of the first chemical of the products in a reaction
 * @param end_bound: the index after the last chemical of the products in a reaction
*/
__global__ void update_chem_amounts_using_products(chem_arr_t chem_arrays, field_t products, int start_bound, int end_bound) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x + start_bound;

    if (i < end_bound) {
        // x normally should be unique to each thread, but in case it isn't guard against race conditions with atomics
        chem_id_t x = products.chem_ids[i];
        atomicAdd(&chem_arrays.chem_amounts[x], products.deltas[i]);
    }
    __syncthreads();
}

/*
 * A kernel function to check of a chemical has triggered a threshold. This is mainly created to avoid
 * a cudaMemcpy call of the chemical amount array for each iteration of the while loop in simulation_master.
 *
 * @param triggered_thresh: an output pointer to one unsigned integer
 * @param within_threshold: an output pointer to one bool
 * @param chem_arrays: a struct containing arrays for the chemical amounts, threshold amounts, and threshold types
 * @param num_chems: the count of all chemical species in the CRN, not the amount of each chemical
*/
__global__ void check_thresholds(bool *within_threshold, chem_arr_t chem_arrays, unsigned int num_chems) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    __shared__ bool block_within_threshold;

    if (threadIdx.x == 0) {
        block_within_threshold = true;
    }

    __syncthreads();
    if (i < num_chems) {
        unsigned int chem_amount = chem_arrays.chem_amounts[i];
        threshold_types thresh_type = chem_arrays.thresh_types[i];
        unsigned int thresh_amount = chem_arrays.thresh_amounts[i];

        switch(thresh_type) {
            // Each thread that detects a triggered threshold will write the same value to block_within_threshold, so race conditions don't matter
            case THRESH_LT:
                if (chem_amount < thresh_amount) {
                    block_within_threshold = false;
                }
                break;
            case THRESH_LE:
                if (chem_amount <= thresh_amount) {
                    block_within_threshold = false;
                }
                break;
            case THRESH_GE:
                if (chem_amount >= thresh_amount) {
                    block_within_threshold = false;
                }
                break;
            case THRESH_GT:
                if (chem_amount > thresh_amount) {
                    block_within_threshold = false;
                }
                break;
            case THRESH_N:
            default:
                break;
        }
    }

    __syncthreads();
    if (threadIdx.x == 0 && !block_within_threshold) {
        // Like before, each block that detects a triggered threshold will write the same value to within_threshold, so race conditions don't matter
        *within_threshold = false;
    }
}

__global__ void count_triggered_thresholds(unsigned int *total_count_triggered_threshs, int *invalid_thresh_type, chem_arr_t chem_arrays, unsigned int num_chems) {
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    __shared__ int local_invalid_thresh_type;

    if (threadIdx.x == 0) {
        local_invalid_thresh_type = 0;
    }

    __syncthreads();
    if (i < num_chems) {
        switch(chem_arrays.thresh_types[i]) {
            case THRESH_LT:
                if (chem_arrays.chem_amounts[i] < chem_arrays.thresh_amounts[i]) {
                    total_count_triggered_threshs[i]++;
                }
                break;
            case THRESH_LE:
                if (chem_arrays.chem_amounts[i] <= chem_arrays.thresh_amounts[i]) {
                    total_count_triggered_threshs[i]++;
                }
            break;
            case THRESH_GE:
                if (chem_arrays.chem_amounts[i] >= chem_arrays.thresh_amounts[i]) {
                    total_count_triggered_threshs[i]++;
                }
                break;
            case THRESH_GT:
                if (chem_arrays.chem_amounts[i] > chem_arrays.thresh_amounts[i]) {
                    total_count_triggered_threshs[i]++;
                }
                break;
            case THRESH_N:
                break;
            default:
                atomicExch(&local_invalid_thresh_type, 1);
        }
    }

    __syncthreads();
    if (threadIdx.x == 0 && local_invalid_thresh_type) {
        atomicExch(invalid_thresh_type, THRESH_CODE_ERR);
    }
}

/*
 * A recursive function that calls scanBlockSumsKernel() and addScanBlocksKernel() to find
 * the inclusive sum of block sums.
 *
 * @param block_sums: an array containing prefix sum of each block
 * @param num_elements: the count of all elements in block_sums
*/
void scanBlockSumArray(float *block_sums, unsigned int num_elements) {
	unsigned int num_blocks = (unsigned int) ceil((double) num_elements / (double) (2*BLOCK_SIZE));

	float *scanned_block_sums;										        // Allocate space for scanned block sums
	cudaMalloc((void**)&scanned_block_sums, num_blocks * sizeof(float));

	dim3 dimGrid = dim3(num_blocks, 1, 1);
    dim3 dimBlock = dim3(BLOCK_SIZE, 1, 1);

	scanBlockSumsKernel<<<dimGrid, dimBlock>>>(scanned_block_sums, block_sums, num_elements);

	if (num_blocks > 1) {													// Recursively call function to scan the scanned block sums
		scanBlockSumArray(scanned_block_sums, num_blocks);
		dim3 newDimBlock = dim3(2 * BLOCK_SIZE, 1, 1);
		addScanBlocksKernel<<<dimGrid, newDimBlock>>>(block_sums, scanned_block_sums, num_elements, num_blocks);
	}

	cudaFree(scanned_block_sums);											// Clean up after yourself
}

/*
 * A function to find the inclusive sum of an arbitrarily-sized input array.
 *
 * @param sum: a pointer to a floating-point variable in the calling function
 * @param out_array: an output array of partial sums
 * @param in_array: an input array of floats to find the prefix sum
 * @param block_sums: the prefix sum of each block
 * @param num_elements: the count of all elements in the output and input arrays
*/
void prescanArray(float *sum, float *out_array, float *in_array, float *block_sums, unsigned int num_elements) {
    //all pointers are to pre-allocated device memory regions

	unsigned int num_blocks = (unsigned int) ceil((double) num_elements / (double) (2*BLOCK_SIZE));

	dim3 dimGrid = dim3(num_blocks, 1, 1);
    dim3 dimBlock = dim3(BLOCK_SIZE, 1, 1);

	prescanArrayKernel<<<dimGrid, dimBlock>>>(out_array, in_array, block_sums, num_elements);

	if (num_blocks <= 1) {													// To make the runtime of one-block array a bit faster
        cudaMemcpy(sum, &out_array[num_elements - 1], sizeof(float), cudaMemcpyDeviceToHost);
		return;
    }

	scanBlockSumArray(block_sums, num_blocks);

	dim3 newDimBlock = dim3(2 * BLOCK_SIZE, 1, 1);
	addScanBlocksKernel<<<dimGrid, newDimBlock>>>(out_array, block_sums, num_elements, num_blocks);
    cudaMemcpy(sum, &out_array[num_elements - 1], sizeof(float), cudaMemcpyDeviceToHost);
}

/*
 * The recursive calling function of kernel function above.
 *
 * @param candidate_reactions: an array of unsigned integers representing a reaction
 * @param num_reactions: the count of all reactions in the CRN
 * @return the chosen reaction
*/
unsigned int choose_reaction(unsigned int *candidate_reactions, unsigned int num_reactions) {
    if (num_reactions <= 1) {
        unsigned int chosen_reaction = 0;
        cudaMemcpy(&chosen_reaction, &candidate_reactions[0], sizeof(unsigned int), cudaMemcpyDeviceToHost);
        return chosen_reaction;
    }

	dim3 dimGrid = dim3(ceil((float) num_reactions / (float) (2*BLOCK_SIZE)), 1, 1);
    dim3 dimBlock = dim3(BLOCK_SIZE, 1, 1);
	choose_reaction_kernel<<<dimGrid, dimBlock>>>(candidate_reactions, num_reactions);

	return choose_reaction(candidate_reactions, num_reactions/2);
}

/*
 * The master function that handles device memory needed for the simulation and transfers to and from it, performs computation needed for
 * said simulation, and some error handling outside of kernel functions.
 *
 * @params: Way too many to list. Just look at the function.
*/
extern "C" simulation_err_t simulation_master(unsigned int *post_trial_chem_amounts_h, unsigned int *total_count_triggered_threshs_h,
                                        crn_t crn_h, sim_params_t sim_params, output_stats_t *out_stats,
                                        unsigned int trial_num) {
    simulation_err_t ret_err = SIMULATION_SUCCESS;

    chem_arr_t chem_arrays_d;
    field_t reactants_d;
    field_t products_d;
    float *rates_d;

    cudaMalloc((void**)&chem_arrays_d.chem_amounts, crn_h.num_chems * sizeof(unsigned int));          // Allocate all device memory needed for simulation.
    cudaMalloc((void**)&chem_arrays_d.thresh_amounts, crn_h.num_chems * sizeof(unsigned int));
    cudaMalloc((void**)&chem_arrays_d.thresh_types, crn_h.num_chems * sizeof(threshold_types));

    cudaMalloc((void**)&reactants_d.chem_ids, crn_h.num_reactants * sizeof(chem_id_t));
    cudaMalloc((void**)&reactants_d.deltas, crn_h.num_reactants * sizeof(unsigned int));
    cudaMalloc((void**)&reactants_d.start_bounds, crn_h.num_reactions * sizeof(unsigned int));
    cudaMalloc((void**)&reactants_d.end_bounds, crn_h.num_reactions * sizeof(unsigned int));

    cudaMalloc((void**)&products_d.chem_ids, crn_h.num_products * sizeof(chem_id_t));
    cudaMalloc((void**)&products_d.deltas, crn_h.num_products * sizeof(unsigned int));
    cudaMalloc((void**)&products_d.start_bounds, crn_h.num_reactions * sizeof(unsigned int));
    cudaMalloc((void**)&products_d.end_bounds, crn_h.num_reactions * sizeof(unsigned int));

    cudaMalloc((void**)&rates_d, crn_h.num_reactions * sizeof(float));

    size_t padded_num_reactions = TILE_SIZE*((crn_h.num_reactions+TILE_SIZE-1)/TILE_SIZE);

    float *propensities;
    float *propensity_inclusive_sums;
    float *propensity_block_sums;

    cudaMalloc((void**)&propensities, padded_num_reactions*sizeof(float));
    cudaMalloc((void**)&propensity_inclusive_sums, padded_num_reactions*sizeof(float));
    cudaMalloc((void**)&propensity_block_sums, crn_h.num_reactions*sizeof(float));

    float *propensity_h = (float*) malloc(crn_h.num_reactions * sizeof(float));

    unsigned int *candidate_reactions;
    cudaMalloc((void**)&candidate_reactions, crn_h.num_reactions * sizeof(unsigned int));

    unsigned int *total_triggered_threshs_d;
    cudaMalloc((void**)&total_triggered_threshs_d, crn_h.num_chems * sizeof(unsigned int));

    bool within_threshold = true;
    bool *within_threshold_d;
    cudaMalloc((void**)&within_threshold_d, sizeof(bool));

    int *is_thresh_type_invalid;
    cudaMalloc((void**)&is_thresh_type_invalid, sizeof(int));

    cudaStream_t chem_arr_stream, reactants_stream, products_stream, rates_stream, propensity_stream, sim_stream;      // Prepare simulation with streams.
    cudaStreamCreate(&chem_arr_stream);
    cudaStreamCreate(&reactants_stream);
    cudaStreamCreate(&products_stream);
    cudaStreamCreate(&rates_stream);
    cudaStreamCreate(&propensity_stream);
    cudaStreamCreate(&sim_stream);

    cudaMemcpyAsync(chem_arrays_d.chem_amounts, crn_h.chem_arrays.chem_amounts, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyHostToDevice, chem_arr_stream);
    cudaMemcpyAsync(chem_arrays_d.thresh_amounts, crn_h.chem_arrays.thresh_amounts, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyHostToDevice, chem_arr_stream);
    cudaMemcpyAsync(chem_arrays_d.thresh_types, crn_h.chem_arrays.thresh_types, crn_h.num_chems * sizeof(threshold_types), cudaMemcpyHostToDevice, chem_arr_stream);

    cudaMemcpyAsync(reactants_d.chem_ids, crn_h.reactants.chem_ids, crn_h.num_reactants * sizeof(chem_id_t), cudaMemcpyHostToDevice, reactants_stream);
    cudaMemcpyAsync(reactants_d.deltas, crn_h.reactants.deltas, crn_h.num_reactants * sizeof(unsigned int), cudaMemcpyHostToDevice, reactants_stream);
    cudaMemcpyAsync(reactants_d.start_bounds, crn_h.reactants.start_bounds, crn_h.num_reactions * sizeof(unsigned int), cudaMemcpyHostToDevice, reactants_stream);
    cudaMemcpyAsync(reactants_d.end_bounds, crn_h.reactants.end_bounds, crn_h.num_reactions * sizeof(unsigned int), cudaMemcpyHostToDevice, reactants_stream);

    cudaMemcpyAsync(products_d.chem_ids, crn_h.products.chem_ids, crn_h.num_products * sizeof(chem_id_t), cudaMemcpyHostToDevice, products_stream);
    cudaMemcpyAsync(products_d.deltas, crn_h.products.deltas, crn_h.num_products * sizeof(unsigned int), cudaMemcpyHostToDevice, products_stream);
    cudaMemcpyAsync(products_d.start_bounds, crn_h.products.start_bounds, crn_h.num_reactions * sizeof(unsigned int), cudaMemcpyHostToDevice, products_stream);
    cudaMemcpyAsync(products_d.end_bounds, crn_h.products.end_bounds, crn_h.num_reactions * sizeof(unsigned int), cudaMemcpyHostToDevice, products_stream);

    cudaMemcpyAsync(rates_d, crn_h.rates, crn_h.num_reactions * sizeof(float), cudaMemcpyHostToDevice, rates_stream);

    cudaMemsetAsync(propensities, 0.0, padded_num_reactions * sizeof(float), propensity_stream);
    cudaMemsetAsync(propensity_inclusive_sums, 0.0, padded_num_reactions * sizeof(float), propensity_stream);
    cudaMemsetAsync(propensity_block_sums, 0.0, crn_h.num_reactions * sizeof(float), propensity_stream);

    cudaMemcpyAsync(total_triggered_threshs_d, total_count_triggered_threshs_h, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyHostToDevice, sim_stream);
    cudaMemsetAsync(is_thresh_type_invalid, 0, sizeof(int), sim_stream);

    cudaDeviceSynchronize();

    if (trial_num <= 0) srand(time(NULL));

    while (within_threshold && (sim_params.max_time <= 0 || out_stats->time_elapsed < sim_params.max_time) && (sim_params.max_steps <= 0 || out_stats->steps_elapsed < sim_params.max_steps)) {
        unsigned int num_blocks = (unsigned int) ceil((float) crn_h.num_reactions/(2*BLOCK_SIZE));

        // Find propensity of each reaction
        dim3 dimGrid = dim3(num_blocks, 1, 1);
        dim3 dimBlock = dim3(2*BLOCK_SIZE, 1, 1);
        calculate_propensity_kernel<<<dimGrid, dimBlock>>>(propensities, chem_arrays_d, reactants_d, rates_d, crn_h.num_reactions);

        // Find propensity inclusive prefix sums
        float propensity_sum = FLT_MAX;
        cudaMemset(propensity_inclusive_sums, 0.0, padded_num_reactions * sizeof(float));
        cudaMemset(propensity_block_sums, 0.0, crn_h.num_reactions * sizeof(float));
        prescanArray(&propensity_sum, propensity_inclusive_sums, propensities, propensity_block_sums, padded_num_reactions);

        if (propensity_sum <= 0.0f) {
            if (sim_params.verbosity_bit_fields & PRINT_TERMINAL) {
                printf("No further reactions are possible\n");
            }
            break;
        } else if (sim_params.verbosity_bit_fields & PRINT_TRACE) {
            cudaMemcpy(propensity_h, propensities, sizeof(float) * crn_h.num_reactions, cudaMemcpyDeviceToHost);

            printf("Prob [");
            for (unsigned j = 0; j < crn_h.num_reactions; j++) {
                printf("%.4f", propensity_h[j] / propensity_sum);
                if (j < crn_h.num_reactions - 1)
                    printf(", ");
            }
            printf("]\n");
        }

        float r = rand()/(float)RAND_MAX;

        // Choose a reaction
        find_reaction_candidates<<<dimGrid, dimBlock>>>(candidate_reactions, propensity_inclusive_sums, propensity_sum, r, crn_h.num_reactions);
        unsigned int chosen_reaction = choose_reaction(candidate_reactions, crn_h.num_reactions);

        if (chosen_reaction >= crn_h.num_reactions) {
            printf("Error: invalid reaction %u has been chosen.\n", chosen_reaction);
            ret_err = INVALID_REACT_CHOSEN;
            break;
        }

        if (sim_params.verbosity_bit_fields & PRINT_TRACE) {
            printf("Chosen reaction: %u\n", chosen_reaction);
        }

        r = rand()/(float)RAND_MAX;

        if (sim_params.max_time >= 0) {
            r = rand()/(double)RAND_MAX;

            double time_step = (-1.0/(double)propensity_sum)*log(r);
            out_stats->time_elapsed += time_step;
            if (sim_params.verbosity_bit_fields & PRINT_TRACE) {
                printf("Chosen time step: %f\n", time_step);
            }
        }

        out_stats->steps_elapsed++;

        unsigned int start_bound_r = crn_h.reactants.start_bounds[chosen_reaction];
        unsigned int end_bound_r = crn_h.reactants.end_bounds[chosen_reaction];

        if (start_bound_r >= end_bound_r) {
            printf("Error: in reactant bound array, start bound is greater than end bound\n");
            ret_err = REAC_BOUND_ARR_ERR;
            break;
        }

        // Remove chemicals that are spent
        dim3 dimGrid_update_r = dim3(ceil((double) (end_bound_r - start_bound_r) / (double) BLOCK_SIZE), 1, 1);
        dim3 dimBlock_update_r = dim3(BLOCK_SIZE, 1, 1);
        update_chem_amounts_using_reactants<<<dimGrid_update_r, dimBlock_update_r>>>(chem_arrays_d, reactants_d, start_bound_r, end_bound_r);

        unsigned int start_bound_p = crn_h.products.start_bounds[chosen_reaction];
        unsigned int end_bound_p = crn_h.products.end_bounds[chosen_reaction];

        if (start_bound_p >= end_bound_p) {
            printf("Error: in product bound array, start bound is greater than end bound\n");
            ret_err = PROD_BOUND_ARR_ERR;
            break;
        }

        // Add chemicals that are produced
        dim3 dimGrid_update_p = dim3(ceil((double) (end_bound_p - start_bound_p) / (double) BLOCK_SIZE), 1, 1);
        dim3 dimBlock_update_p = dim3(BLOCK_SIZE, 1, 1);
        update_chem_amounts_using_products<<<dimGrid_update_p, dimBlock_update_p>>>(chem_arrays_d, products_d, start_bound_p, end_bound_p);

        if (sim_params.verbosity_bit_fields & PRINT_STATES) {
            cudaMemcpy(post_trial_chem_amounts_h, chem_arrays_d.chem_amounts, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyDeviceToHost);
            printf("S [");
            for (unsigned j = 0; j < crn_h.num_chems; j++) {
                printf("%3u", post_trial_chem_amounts_h[j]);
               if (j < crn_h.num_chems - 1) printf(", ");
            }
            printf("]\n");
        }

        // Check if a threshold has been triggered
        dim3 dimGrid_thresh = dim3(ceil((double) crn_h.num_chems / (double) BLOCK_SIZE), 1, 1);
        dim3 dimBlock_thresh = dim3(BLOCK_SIZE, 1, 1);
        cudaMemset(within_threshold_d, true, sizeof(bool));
        check_thresholds<<<dimGrid_thresh, dimBlock_thresh>>>(within_threshold_d, chem_arrays_d, crn_h.num_chems);

        cudaMemcpy(&within_threshold, within_threshold_d, sizeof(bool), cudaMemcpyDeviceToHost);

        if (!within_threshold) {
            count_triggered_thresholds<<<dimGrid_thresh, dimBlock_thresh>>>(total_triggered_threshs_d, is_thresh_type_invalid, chem_arrays_d, crn_h.num_chems);
            cudaMemcpy(&ret_err, is_thresh_type_invalid, sizeof(int), cudaMemcpyDeviceToHost);

            if (ret_err) {
                printf("Error: invalid threshold code. Please check your .in file for errors.");
                break;
            }
        }
    }

    // If print states was enabled, then an up-to-date host-side copy of the chem amounts already exists. Therefore, don't copy it again.
    if (ret_err == SIMULATION_SUCCESS) {
        cudaMemcpy(post_trial_chem_amounts_h, chem_arrays_d.chem_amounts, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyDeviceToHost);

        if (!within_threshold) {
            cudaMemcpy(total_count_triggered_threshs_h, total_triggered_threshs_d, crn_h.num_chems * sizeof(unsigned int), cudaMemcpyDeviceToHost);

            if (sim_params.verbosity_bit_fields & PRINT_STATES) {
                printf("State after threshold [");
                for (unsigned j = 0; j < crn_h.num_chems; j++) {
                    printf("%3u", post_trial_chem_amounts_h[j]);
                    if (j < crn_h.num_chems - 1) printf(", ");
                }
                printf("]\n");
            }
            if (sim_params.verbosity_bit_fields & PRINT_TERMINAL) {
                printf("State exceeds threshold.\n");
            }
        }

        if (out_stats->steps_elapsed >= sim_params.max_steps && sim_params.max_steps > 0) {
            printf("Max steps has been reached.\n");
        }
    }

    cudaStreamDestroy(chem_arr_stream);
    cudaStreamDestroy(reactants_stream);
    cudaStreamDestroy(products_stream);
    cudaStreamDestroy(rates_stream);
    cudaStreamDestroy(propensity_stream);
    cudaStreamDestroy(sim_stream);

    free(propensity_h);

    cudaFree(within_threshold_d);
    cudaFree(propensities);
    cudaFree(propensity_inclusive_sums);
    cudaFree(propensity_block_sums);

    cudaFree(candidate_reactions);

    cudaFree(chem_arrays_d.chem_amounts);
    cudaFree(chem_arrays_d.thresh_amounts);
    cudaFree(chem_arrays_d.thresh_types);

    cudaFree(reactants_d.chem_ids);
    cudaFree(reactants_d.deltas);
    cudaFree(reactants_d.start_bounds);
    cudaFree(reactants_d.end_bounds);

    cudaFree(products_d.chem_ids);
    cudaFree(products_d.deltas);
    cudaFree(products_d.start_bounds);
    cudaFree(products_d.end_bounds);

    cudaFree(rates_d);

    return ret_err;
}
