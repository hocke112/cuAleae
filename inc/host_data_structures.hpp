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
 * This file contains definitions of structs used to hold host-side data structures, improving upon the ones in aleae.h of Aleae's source code.
*/

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

    float rate = 0;
} typedef reaction_t;

struct CRN {
    std::vector<chem_t> chems;
    std::vector<std::unique_ptr<reaction_t>> reactions;
    std::vector<thresh_t> thresholds;
};

#endif