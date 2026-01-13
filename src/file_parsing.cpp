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
 * This file contains the logic to parse the .in and .r input files.
*/

#include <iostream>
#include <fstream>
#include <regex>

#include "file_parsing.hpp"

/*
 * A function to parse a .in input file containing the initial state of the CRN.
 *
 * @param in_filename: a string containing the filename
 * @param chem_str_to_id: a mapping of chemical names to their chem IDs
 * @param crn: the chemical reaction network
 * @return: 1 of an error has occurred or 0 otherwise
*/
int parse_in_input_file(const std::string &in_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn) {
    std::ifstream in_file(in_filename);
    if (in_filename.substr(in_filename.length() - 3, 3) != ".in") {
        std::cerr << "Error: Input .in file type is invalid for .in file." << std::endl;
        return 1;
    }

    if (!in_file.is_open()) {
        std::cerr << "Error: Input .in file failed to be opened" << std::endl;
        return 1;
    }

    std::regex init_line_re("([A-Za-z'][A-Za-z0-9.'_]*)[ ]+([0-9]+)[ ]+(LE|LT|GE|GT|N)([ ]+([0-9]*))*");
    std::regex cr_char_re("\\r");

    std::string temp("");
    std::string chem_name("");
    std::string chem_amount_str("");
    std::string thresh_type_str("");
    std::string thresh_amount_str("");

    unsigned int lineno = 0;
    while (getline(in_file, temp)) {
        temp = std::regex_replace(temp, cr_char_re, "");

        chem_name = std::regex_replace(temp, init_line_re, "$1");
        chem_amount_str = std::regex_replace(temp, init_line_re, "$2");
        thresh_type_str = std::regex_replace(temp, init_line_re, "$3");
        thresh_amount_str = std::regex_replace(temp, init_line_re, "$5");

        if (chem_name == temp || chem_amount_str == temp || thresh_type_str == temp) {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Syntax error at line " << lineno + 1 << std::endl;
            return 1;
        } else if (thresh_amount_str == "") {
            thresh_amount_str = "0";
        }

        unsigned int chem_amount = (unsigned int) std::stoul(chem_amount_str);

        threshold_types thresh_type = THRESH_N;

        if (thresh_type_str == "LT") {
            thresh_type = THRESH_LT;
        } else if (thresh_type_str == "LE") {
            thresh_type = THRESH_LE;
        } else if (thresh_type_str == "GE") {
            thresh_type = THRESH_GE;
        } else if (thresh_type_str == "GT") {
            thresh_type = THRESH_GT;
        }

        unsigned int thresh_amount = (unsigned int) std::stoul(thresh_amount_str);

        if (thresh_type == THRESH_N && thresh_amount > 0 || thresh_type < THRESH_N && thresh_amount == 0) {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Syntax error at line " << lineno + 1 << std::endl;
            return 1;
        }

        if (chem_str_to_id.count(chem_name) > 0) {
            std::cerr << "In " << in_filename << std::endl;
            std::cerr << "Duplicate chemical at line " << lineno + 1 << std::endl;
            return 1;
        } else {
            chem_str_to_id.insert(std::pair(chem_name, lineno));
        }

        crn.chems.emplace_back(chem_t(chem_name, chem_amount));
        crn.thresholds.emplace_back(thresh_t(thresh_type, thresh_amount));

        lineno++;
    }

    in_file.close();
    return 0;
}

/*
 * A function to parse a .r input file containing the reactions in a CRN.
 *
 * @param r_filename: a string containing the filename
 * @param chem_str_to_id: a mapping of chemical names to their chem IDs
 * @param crn: the chemical reaction network
 * @return: 1 of an error has occurred or 0 otherwise
*/
int parse_r_input_file(const std::string &r_filename, const std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn) {
    std::ifstream r_file(r_filename);
    if (r_filename.substr(r_filename.length() - 2, 2) != ".r") {
        std::cerr << "Error: Input .r file type is invalid for .in file." << std::endl;
        return 1;
    }

    if (!r_file.is_open()) {
        std::cerr << "Error: Input .r file failed to be opened" << std::endl;
        return 1;
    }

    std::regex integer_re("([1-9][0-9]*)");
    std::regex equation_re("([^:]+)[ ]*:[ ]*([^:]+)[ ]*:[ ]*([0-9]+([.][0-9]+((E|e)[-]{0,1}[0-9]+)*)*)[ \\r]*$");
    std::regex term_re("([A-Za-z'][A-Za-z0-9.'_]*) ([1-9][0-9]*)");

    std::regex cr_char_re("\\r");

    unsigned int reaction_no = 0;
    std::string reactant_str("");
    std::string product_str("");
    std::string rate_str("");
    std::string temp("");
    while (getline(r_file, temp)) {
        temp = std::regex_replace(temp, cr_char_re, "");

        reactant_str = std::regex_replace(temp, equation_re, "$1");
        product_str = std::regex_replace(temp, equation_re, "$2");
        rate_str = std::regex_replace(temp, equation_re, "$3");

        if (reactant_str == temp || product_str == temp || rate_str == temp) {
            std::cerr << "In " << r_filename << std::endl;
            std::cerr << "Syntax error at line " << reaction_no + 1 << std::endl;
            return 1;
        }

        if (rate_str == "0") {
            std::cerr << "In " << r_filename << std::endl;
            std::cerr << "Error at line " << reaction_no + 1 << std::endl;
            std::cerr << "Reaction rate must be a positive non-zero value." << std::endl;
            return 1;
        }

        crn.reactions.emplace_back(std::make_unique<reaction_t>());
        crn.reactions[reaction_no]->rate = std::stof(rate_str);

        for (std::sregex_iterator it_reactant = std::sregex_iterator(reactant_str.begin(), reactant_str.end(), term_re);
                it_reactant != std::sregex_iterator(); ++it_reactant){
            std::string cur_term = std::smatch(*it_reactant).str();

            std::string chem_str = std::regex_replace(cur_term, term_re, "$1");
            std::string coeff_st = std::regex_replace(cur_term, term_re, "$2");

            if (!chem_str_to_id.count(chem_str)) {
                std::cerr << "Missing chemical at line " << reaction_no + 1 << std::endl;
                return 1;
            }

            crn.reactions[reaction_no]->reactants.emplace_back(chem_str_to_id.at(chem_str));
            crn.reactions[reaction_no]->reactant_deltas.emplace_back((unsigned int) std::stoul(coeff_st));
        }

        for (std::sregex_iterator it_product = std::sregex_iterator(product_str.begin(), product_str.end(), term_re);
                it_product != std::sregex_iterator(); ++it_product){
            std::string cur_term = std::smatch(*it_product).str();

            std::string chem_str = std::regex_replace(cur_term, term_re, "$1");
            std::string coeff_st = std::regex_replace(cur_term, term_re, "$2");

            if (!chem_str_to_id.count(chem_str)) {
                std::cerr << "In " << r_filename << std::endl;
                std::cerr << "Missing chemical at line " << reaction_no + 1 << std::endl;
                return 1;
            }

            crn.reactions[reaction_no]->products.emplace_back(chem_str_to_id.at(chem_str));
            crn.reactions[reaction_no]->product_deltas.emplace_back((unsigned int) std::stoul(coeff_st));
        }

        reaction_no++;
    }

    r_file.close();
    return 0;
}