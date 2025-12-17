#ifndef __HOST_DATA_STRUCTURES_HPP
#define __HOST_DATA_STRUCTURES_HPP

#include <string>
#include <vector>
#include <memory>
#include "misc_def.h"

struct chemical {
    std::string name;
    unsigned int amount = 0;

    chemical(std::string name, unsigned int amount) {
        this->name = name;
        this->amount = amount;
    }

    chemical() {
        this->name = "";
    }
} typedef chem_t;

struct threshold {
    threshold_types type = THRESH_N;
    unsigned int amount = 0;

    threshold(threshold_types type, unsigned int amount) {
        this->type = type;
        this->amount = amount;
    }

    threshold() {
        return;
    }
} typedef thresh_t;

struct reaction {
    std::vector<chem_id_t> reactants;
    std::vector<unsigned int> reactant_deltas;      // Coefficients in a chemical reaction

    std::vector<chem_id_t> products;
    std::vector<unsigned int> product_deltas;       // Coefficients in a chemical reaction

    double rate = 0;
} typedef reaction_t;

struct CRN {
    std::vector<chem_t> chems;
    std::vector<std::unique_ptr<reaction_t>> reactions;
    std::vector<thresh_t> thresholds;
};

#endif