#include <iostream>
#include <regex>
#include <chrono>
#include <iomanip>
#include <unordered_map>
#include "host_data_structures.hpp"
#include "device_data_structures.h"
#include "file_parsing.hpp"
#include "simulation.h"

/*
 * Prints all the contents of each of the dynamically-allocated arrays in the CRN struct.
 *
 * @param: a struct representing the CRN.
*/
int print_crn_contents(const CRN &crn) {
    std::cout << "Initial Quantities and Thresholds\n";
    for (unsigned int i = 0; i < crn.chems.size(); ++i) {
        std::cout << crn.chems[i].name << " " << crn.chems[i].amount;
        std::cout << " ";

        threshold_types t = crn.thresholds[i].type;
        switch (t) {
            case THRESH_LT:
                std::cout << "(< ";
                break;
            case THRESH_LE:
                std::cout << "(<= ";
                break;
            case THRESH_GE:
                std::cout << "(>= ";
                break;
            case THRESH_GT:
                std::cout << "(> ";
                break;
            case THRESH_N:
                break;
            default:
                std::cerr << "error: invalid threshold code" << std::endl;
                return 1;
        }

        std::cout << crn.thresholds[i].amount << ")";
        if (i  < crn.chems.size() - 1) std::cout << ", ";
    }

    std::cout << "\n\n";
    std::cout << "Reactions\n";
    for (unsigned int k = 0; k < crn.reactions.size(); ++k) {
        for (unsigned int l = 0; l < crn.reactions[k]->reactants.size(); ++l) {
            std::cout << crn.chems[crn.reactions[k]->reactants[l]].name << " " << crn.reactions[k]->reactant_deltas[l] << " ";
        }
        std::cout << ": ";
        for (unsigned int m = 0; m < crn.reactions[k]->reactants.size(); ++m) {
            std::cout << crn.chems[crn.reactions[k]->reactants[m]].name << " " << crn.reactions[k]->reactant_deltas[m] << " ";
        }
        std::cout << ": " << crn.reactions[k]->rate << "\n";
    }
    std::cout << "\n";

    return 0;
}

int main (int argc, char *argv[]) {
    if (argc != 7) {
        std::cerr << "Missing arguments." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
     }

    std::regex integer_re("([-]{0,1}[0-9]+)");
    std::smatch sm;
    std::string trial_str = argv[3];
    if (!std::regex_search(trial_str, sm, integer_re)) {
        std::cerr << "Value of trial must be a positive non-zero integer." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    std::string time_str = argv[4];
    if (!std::regex_search(time_str, sm, integer_re)) {
        std::cerr << "Time value must be an integer." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    std::string verbosity_str = argv[5];
    if (!std::regex_search(verbosity_str, sm, integer_re)) {
        std::cerr << "Verbosity value must be an integer." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    std::string steps_str = argv[6];
    if (!std::regex_search(steps_str, sm, integer_re)) {
        std::cerr << "Steps value must be an integer." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    std::string in_filename(argv[1]);
    std::string r_filename(argv[2]);

    int num_trials = std::stoi(argv[3]);
    int num_time = std::stoi(argv[4]);
    int verbosity = std::stoi(argv[5]);
    int num_steps = std::stoi(argv[6]);

    if (num_trials <= 0) {
        std::cerr << "Value of trial must be an a positive non-zero value." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    if (num_time < -1) {
        std::cerr << "Verbosity value must be greater than or equal to -1." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    if (verbosity > 15 || verbosity < 0) {
        std::cerr << "Verbosity value must be between 0 and 15 inclusive." << std::endl;
        std::cerr << "usage: " << argv[0] << " <file1: state> <file2: reactions> <trials> <time> <verbosity> <steps>" << std::endl;
        return 1;
    }

    std::map<std::string, chem_id_t> chem_str_to_id;

    //Intialize Chemical Reaction Network (crn) read in and parsed by input files
    CRN crn;

    if (parse_in_input_file(in_filename, chem_str_to_id, crn)) {
        return 1;
    }

    if (parse_r_input_file(r_filename, chem_str_to_id, crn)){
        return 1;
    }

    const auto start_prog = std::chrono::high_resolution_clock::now();

    // define output_stats_t
    output_stats_t out_stats;

    // define sim_params_t
    sim_params_t sim_params;
    sim_params.verbosity_bit_fields = verbosity;
    sim_params.max_time = num_time;


    //total number of chemicals in network
    unsigned int num_chems = crn.chems.size();
    //total number of reactions in network
    unsigned int num_reactions = crn.reactions.size();


    // This counts all of the reactants an totals products in all reactions
    size_t total_reactants = 0;
    size_t total_products  = 0;

    //for loops loops through all reactions and adds these amounts to our total variables
    for (auto& rptr : crn.reactions) {
        total_reactants += rptr->reactants.size();
        total_products  += rptr->products.size();
    }

    //Allocates all reactants
    field_t reactants_h;
    reactants_h.chem_ids     = (chem_id_t*)malloc(total_reactants * sizeof(chem_id_t));
    reactants_h.deltas       = (unsigned int*)malloc(total_reactants * sizeof(unsigned int));
    reactants_h.start_bounds = (unsigned int*)malloc(num_reactions * sizeof(unsigned int));
    reactants_h.end_bounds   = (unsigned int*)malloc(num_reactions * sizeof(unsigned int));

    // Allocates all products
    field_t products_h;
    products_h.chem_ids     = (chem_id_t*)malloc(total_products * sizeof(chem_id_t));
    products_h.deltas       = (unsigned int*)malloc(total_products * sizeof(unsigned int));
    products_h.start_bounds = (unsigned int*)malloc(num_reactions * sizeof(unsigned int));
    products_h.end_bounds   = (unsigned int*)malloc(num_reactions * sizeof(unsigned int));

    //populate the reactant host struct and product host struct
    size_t r_idx = 0; //reaction index
    size_t p_idx = 0; // product index
    for (unsigned i = 0; i < num_reactions; i++) {
        const reaction_t* r = crn.reactions[i].get();

        if (r->reactant_deltas.size() != r->reactants.size()) {
            std::cerr << "Reactants are not paired with coefficients properly. Check " << argv[2] << " for any syntax errors." << std::endl;
            return -1;
        } else if (r->reactant_deltas.size() == 0 || r->reactants.size() == 0) {
            std::cerr << "There's a syntax error that was not detected earlier. Check " << argv[2] << " for any syntax errors." << std::endl;
            return -1;
        }

        // fiiling reactants
        reactants_h.start_bounds[i] = r_idx;
        for (size_t j = 0; j < r->reactants.size(); j++) {
            reactants_h.chem_ids[r_idx] = r->reactants[j];
            reactants_h.deltas[r_idx]   = r->reactant_deltas[j];
            r_idx++;
        }
        reactants_h.end_bounds[i] = r_idx;

        if (r->product_deltas.size() != r->products.size()) {
            std::cerr << "Products are not paired with coefficients properly. Check " << argv[2] << " for any syntax errors." << std::endl;
            return -1;
        } else if (r->product_deltas.size() == 0 || r->products.size() == 0) {
            std::cerr << "There's a syntax error that was not detected earlier. Check " << argv[2] << " for any syntax errors." << std::endl;
            return -1;
        }

        // filling products
        products_h.start_bounds[i] = p_idx;
        for (size_t j = 0; j < r->products.size(); j++) {
            products_h.chem_ids[p_idx] = r->products[j];
            products_h.deltas[p_idx]   = r->product_deltas[j];
            p_idx++;
        }
        products_h.end_bounds[i] = p_idx;
    }

    //Put all the data in the struct to hold all of the arrays of chemicals/substances in the network
    chem_arr_t chem_arrays_h;
    unsigned int *post_trial_chem_amounts;

    //allocates components of the chem array for chem amounts and threshhold information
    chem_arrays_h.chem_amounts   = (unsigned int*)malloc(num_chems * sizeof(unsigned int));
    chem_arrays_h.thresh_amounts = (unsigned int*)malloc(num_chems * sizeof(unsigned int));
    chem_arrays_h.thresh_types   = (threshold_types*)malloc(num_chems * sizeof(threshold_types));

    post_trial_chem_amounts      = (unsigned int*)malloc(num_chems * sizeof(unsigned int));

    if (num_chems != crn.thresholds.size()) {
        std::cerr << "Insufficient amount of thresholds detected. Check " << argv[1] << " for any missing thresholds." << std::endl;
        return -1;
    }

    for (unsigned i = 0; i < num_chems; i++) {
        chem_arrays_h.chem_amounts[i]   = crn.chems[i].amount;
        chem_arrays_h.thresh_amounts[i] = crn.thresholds[i].amount;
        chem_arrays_h.thresh_types[i]   = crn.thresholds[i].type;
    }

    //allocates memory based on the total number of rections to have a collection of all the rates
    float* rates_h = (float*)malloc(num_reactions * sizeof(float));
    for (unsigned i = 0; i < num_reactions; i++)
        rates_h[i] = static_cast<float>(crn.reactions[i]->rate);

    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(4);
    if (print_crn_contents(crn)) {
        return -1;
    }

    std::vector<double> avg_post_trial_chem_amounts;
    std::vector<unsigned int> total_triggered_threshs;
    avg_post_trial_chem_amounts.resize(num_chems);
    total_triggered_threshs.resize(num_chems);
    chem_id_t triggered_thresh = UINT32_MAX;


    std::vector<std::unordered_map<unsigned int, int>> z;
    std::vector<std::unordered_map<unsigned int, double>> Prob;

    // std::vector<std::vector<int>> z;
    // std::vector<std::vector<double>> Prob;

    z.resize(num_chems);
    // for (unsigned int i = 0; i < num_chems; ++i) {
    //     z[i].resize(1000);
    // }
    Prob.resize(num_chems);
    // for (unsigned int i = 0; i < num_chems; ++i) {
    //     Prob[i].resize(1000);
    // }

    int err = 0;
    output_stats_t trial_stats;
    bool within_threshold = true;
    for(int i = 0; i < num_trials; ++i){
        //Calls the simulation master for the given trial i
        if (sim_params.verbosity_bit_fields & PRINT_TRIALS) {
            std::cout << "Trial " << i << "\n";
        }

        out_stats.steps_elapsed = 0;
        out_stats.time_elapsed = 0;

        const auto start_trial = std::chrono::high_resolution_clock::now();
        simulation_master(post_trial_chem_amounts, chem_arrays_h, reactants_h, products_h, rates_h,
                                            num_chems, num_reactions, total_reactants, total_products,
                                            sim_params, &out_stats, num_steps, &triggered_thresh, &within_threshold, i, &err);

        if (err) {
            return -1;
        }

        trial_stats.steps_elapsed += out_stats.steps_elapsed;
        trial_stats.time_elapsed += out_stats.time_elapsed;

        for (int j = 0; j < num_chems; j++) {
            avg_post_trial_chem_amounts[j] += (double) post_trial_chem_amounts[j] / (double) num_trials;
            if (!z[j].count(post_trial_chem_amounts[j])) z[j][post_trial_chem_amounts[j]] = 0;
            z[j][post_trial_chem_amounts[j]]++;
        }

        for (unsigned i = 0; i < num_chems && !within_threshold; ++i) {
            switch(crn.thresholds[i].type) {
                case THRESH_LT:
                    if (post_trial_chem_amounts[i] < crn.thresholds[i].amount) {
                        if (sim_params.verbosity_bit_fields & PRINT_TRIALS)
                            std::cout << crn.chems[i].name << " < " << crn.thresholds[i].amount << "\n";
                        total_triggered_threshs[i]++;
                    }
                    break;
                case THRESH_LE:
                    if (post_trial_chem_amounts[i] <= crn.thresholds[i].amount) {
                        if (sim_params.verbosity_bit_fields & PRINT_TRIALS)
                            std::cout << crn.chems[i].name << " <= " << crn.thresholds[i].amount << "\n";
                        total_triggered_threshs[i]++;
                    }
                break;
                case THRESH_GE:
                    if (post_trial_chem_amounts[i] >= crn.thresholds[i].amount) {
                        if (sim_params.verbosity_bit_fields & PRINT_TRIALS)
                            std::cout << crn.chems[i].name << " >= " << crn.thresholds[i].amount << "\n";
                        total_triggered_threshs[i]++;
                    }
                    break;
                case THRESH_GT:
                    if (post_trial_chem_amounts[i] > crn.thresholds[i].amount) {
                        if (sim_params.verbosity_bit_fields & PRINT_TRIALS)
                            std::cout << crn.chems[i].name << " > " << crn.thresholds[i].amount << "\n";
                        total_triggered_threshs[i]++;
                    }
                    break;
                case THRESH_N:
                    break;
                default:
                    std::cerr << "Error: invalid threshold code " << within_threshold << std::endl;
                    free(chem_arrays_h.chem_amounts);
                    free(chem_arrays_h.thresh_amounts);
                    free(chem_arrays_h.thresh_types);
                    free(reactants_h.chem_ids);
                    free(reactants_h.deltas);
                    free(reactants_h.start_bounds);
                    free(reactants_h.end_bounds);
                    free(products_h.chem_ids);
                    free(products_h.deltas);
                    free(products_h.start_bounds);
                    free(products_h.end_bounds);
                    free(rates_h);
                    return -1;
                }
        }

        if (sim_params.verbosity_bit_fields & PRINT_TRIALS) {
            std::cout << "\nTrial stats:\n";

            if (!within_threshold && triggered_thresh >= num_chems) {
                std::cerr << "Error: Invalid threshold found" << std::endl;
                free(chem_arrays_h.chem_amounts);
                free(chem_arrays_h.thresh_amounts);
                free(chem_arrays_h.thresh_types);
                free(reactants_h.chem_ids);
                free(reactants_h.deltas);
                free(reactants_h.start_bounds);
                free(reactants_h.end_bounds);
                free(products_h.chem_ids);
                free(products_h.deltas);
                free(products_h.start_bounds);
                free(products_h.end_bounds);
                free(rates_h);
                return -1;
            }

            std::cout << "Events  " << out_stats.steps_elapsed << "\n";
            if (sim_params.max_time >= 0) {
               std::cout << "Time    " << out_stats.time_elapsed  << "\n";
            }

            const auto end_trial = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double, std::ratio<1,1>> trial_duration = end_trial - start_trial;

            std::cout << "Runtime " << trial_duration.count() << " s\n";
            if (i < num_trials - 1) std::cout << "\n";
        }
    }

    std::vector<double> sum;
    std::vector<std::unordered_map<unsigned int, double>> diff;
    // std::vector<std::vector<double>> diff;
    std::vector<double> var;
    std::vector<double> count;

    sum.resize(num_chems);
    diff.resize(num_chems);
    // for (unsigned int i = 0; i < num_chems; ++i) {
    //     diff[i].resize(1000);
    // }
    var.resize(num_chems);
    count.resize(num_chems);

    // Compute probabilities
	for(unsigned int i = 0; i < num_chems; ++i) {
        for (const std::pair<unsigned int, int> z_col : z[i]) {
            Prob[i][z_col.first] = z_col.second / (double) num_trials;

            if (z_col.second != 0) {
                sum[i] += z_col.first * Prob[i][z_col.first];
                count[i]+=1;
            }
            if (Prob[i][z_col.first] != 0) {
                std::cout << " Probability of (" << crn.chems[i].name << " = "<< z_col.first << ") = "<< Prob[i][z_col.first] << "\n";
            }
        }
		// for(unsigned int j = 0; j<1000; ++j) {
		// 	Prob[i][j] = z[i][j] / (double) num_trials;
		// 	if(z[i][j] != 0) {
		// 		sum[i]+=j*Prob[i][j];
		// 		count[i]+=1;
		// 	}
		// 	if (Prob[i][j] != 0) {
		// 		std::cout <<" Probability of ("<< crn.chems[i].name << " = "<< j << ") = "<< Prob[i][j] << "\n";
		// 	}
		// }
	}

    for(int i = 0; i < num_chems; ++i) {
		std::cout<< " Mean of probability distribution of " << crn.chems[i].name << " = "<< sum[i] << "\n";
	}

	for(int i = 0; i < num_chems; ++i) {
        for (const std::pair<unsigned int, double> d_col : Prob[i]) {
            diff[i][d_col.first] = (d_col.first - sum[i]) * (d_col.first - sum[i]);
            diff[i][d_col.first] *= Prob[i][d_col.first];
            var[i] += diff[i][d_col.first];
        }
		// for(int j=0;j<1000;++j) {
		// 	diff[i][j]=(j-sum[i]) * (j-sum[i]);
		// 	diff[i][j]=diff[i][j] * Prob[i][j];
		// 	var[i]=var[i]+diff[i][j];
		// }
	}

	for(int i = 0; i < num_chems; ++i) {
		std::cout << " Variance of probability distribution of "<< crn.chems[i].name << " = "<< var[i] << "\n";
	}

    std::cout << "\nSimulation stats:\n";
    std::cout << "avg ";
    std::cout << "[";
    int i = 0;
    for (double chem_amounts: avg_post_trial_chem_amounts) {
        std::cout << crn.chems[i].name << ": " << avg_post_trial_chem_amounts[i];
        if (i < avg_post_trial_chem_amounts.size() - 1)
            std::cout << ", ";
        ++i;
    }
    std::cout << "]" << std::endl;

    for (unsigned i = 0; i < num_chems; i++) {
        if (crn.thresholds[i].type == THRESH_N) {
            continue;
        } else {
            std::cout << crn.chems[i].name;
        }
        switch(crn.thresholds[i].type) {
            case THRESH_LT:
                std::cout << " <  ";
                break;
            case THRESH_LE:
                std::cout << " <= ";
                break;
            case THRESH_GE:
                std::cout << " >= ";
                break;
            case THRESH_GT:
                std::cout << " >  ";
                break;
            case THRESH_N:
                break;
            default:
                std::cerr << "Error: invalid threshold code" << std::endl;
                free(chem_arrays_h.chem_amounts);
                free(chem_arrays_h.thresh_amounts);
                free(chem_arrays_h.thresh_types);
                free(reactants_h.chem_ids);
                free(reactants_h.deltas);
                free(reactants_h.start_bounds);
                free(reactants_h.end_bounds);
                free(products_h.chem_ids);
                free(products_h.deltas);
                free(products_h.start_bounds);
                free(products_h.end_bounds);
                free(rates_h);
                return 1;
            }
        std::cout << crn.thresholds[i].amount << ": " << total_triggered_threshs[i] << " (";
        std::cout << ((double) total_triggered_threshs[i]/(double)num_trials)*100 << "%)\n";
    }

    const auto end_prog = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::ratio<1,1>> trial_duration = end_prog - start_prog;
    std::cout << "avg   events/trial " << trial_stats.steps_elapsed/(double)num_trials << "\n";
    if (sim_params.max_time >= 0) {
        std::cout << "avg   time/trial   " << trial_stats.time_elapsed / (double)num_trials << "\n";
    }
    std::cout << "avg   runtime      " << trial_duration.count() / (double)num_trials << " s\n";
    std::cout << "total runtime      " << trial_duration.count() << " s\n";

    //frees at the end and deletes
    free(chem_arrays_h.chem_amounts);
    free(chem_arrays_h.thresh_amounts);
    free(chem_arrays_h.thresh_types);
    free(reactants_h.chem_ids);
    free(reactants_h.deltas);
    free(reactants_h.start_bounds);
    free(reactants_h.end_bounds);
    free(products_h.chem_ids);
    free(products_h.deltas);
    free(products_h.start_bounds);
    free(products_h.end_bounds);
    free(rates_h);

    return 0;
}