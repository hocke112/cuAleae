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

#define REACTANT_FIELD 0
#define PRODUCT_FIELD 1
#define RATE_FIELD 2

/*
 * A helper function to remove whitespace at beginning of a string.
 * @param s: a string to be trimmed
*/
void trim_left(std::string &s) {
    if (s.length() < 1)
        return;

    unsigned int i = 0;
    while (i <= s.length() && s[i] == ' ') {
        ++i;
    }

    s.erase(0, i);
}

/*
 * A helper function to remove whitespace at end of a string.
 * @param s: a string to be trimmed
*/
void trim_right(std::string &s) {
    if (s.length() < 1)
        return;

    unsigned int i = s.length() - 1;
    while (i > 0 && s[i] == ' ') {
        --i;
    }

    if (s[i] != ' ')
        ++i;

    s.erase(i, s.length());
}

/*
 * A helper function to split a string into a list of tokens.
 *
 * @param tokens: an empty list of tokens that will be filled
 * @param s: an input string to get tokens from
 * @char delim: a delimiting character to split the string
 * @return num_split: the number of splits
*/
unsigned int split_by(std::vector<std::string> &tokens, const std::string &s, char delim) {
    unsigned int num_splits = 0;

    unsigned int last_idx = 0;
    for (int i = 0; i < s.length(); ++i) {
        if (s[i] == delim) {
            tokens.emplace_back(s.substr(last_idx, i - last_idx));
            last_idx = i + 1;
            num_splits++;
        }
    }

    tokens.emplace_back(s.substr(last_idx, s.length() - last_idx));

    return num_splits;
}

void print_error_msg_at_line(const std::string filename, unsigned int line_no, const std::string msg) {
    std::cerr << "In " << filename << std::endl;
    std::cerr << "Error at line " << line_no + 1 << std::endl;
    std::cerr << msg << std::endl;
}

/*
 * A function to parse a .in input file containing the initial state of the CRN.
 *
 * @param in_filename: a string containing the filename
 * @param chem_str_to_id: a mapping of chemical names to their chem IDs
 * @param crn: the chemical reaction network
 * @return: 1 of an error has occurred or 0 otherwise
*/
int parse_in_input_file(const std::string &in_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn) {
    if (in_filename.length() < 3) {
        std::cerr << "Error: Input file name " << in_filename << " is invalid." << std::endl;
        return 1;
    } else if (in_filename.substr(in_filename.length() - 3, 3) != ".in") {
        std::cerr << "Error: Input file type of " << in_filename << " is invalid." << std::endl;
        return 1;
    }

    std::ifstream in_file(in_filename);
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
 * A function that fills the reactant or product of a reaction within a CRN.
 *
 * @param reaction_ptr: a unique pointer to a reaction
 * @param chem_str_to_id: a mapping of chemical names to their chem IDs
 * @param token: a token containing the reactant or product field
 * @param field: a field represented as a number
 * @param r_filename: a string containing the filename
 * @param reaction_no: an integer corresponding to a line in the .r file
 * @return 1 of an error has occurred or 0 otherwise
*/
int fill_reaction_field(std::unique_ptr<reaction_t> &reaction_ptr, const std::map<std::string, chem_id_t> &chem_str_to_id, const std::string &token,
                            unsigned int field, const std::string &r_filename, unsigned int reaction_no) {
    static std::regex term_re("([A-Za-z'][A-Za-z0-9.'_]*) ([1-9][0-9]*)");
    static std::regex field_re("^([^: ]+ [0-9]+[ ]*)+$");

    std::smatch sm;
    if (!std::regex_match(token, sm, field_re)) {
        if (field) {
            print_error_msg_at_line(r_filename, reaction_no, "Term in product side of reaction is formatted incorrectly.");
        } else {
            print_error_msg_at_line(r_filename, reaction_no, "Term in reactant side of reaction is formatted incorrectly.");
        }
        return 1;
    }

    for (std::sregex_iterator it_field = std::sregex_iterator(token.begin(), token.end(), term_re);
            it_field != std::sregex_iterator(); ++it_field){
        std::string cur_term = std::smatch(*it_field).str();

        std::string chem_str = std::regex_replace(cur_term, term_re, "$1");
        std::string coeff_st = std::regex_replace(cur_term, term_re, "$2");

        if (!chem_str_to_id.count(chem_str)) {
            print_error_msg_at_line(r_filename, reaction_no, chem_str + " is missing in .in file.");
            return 1;
        }

        if (field) {
            reaction_ptr->products.emplace_back(chem_str_to_id.at(chem_str));
            reaction_ptr->product_deltas.emplace_back((unsigned int) std::stoul(coeff_st));
        } else {
            reaction_ptr->reactants.emplace_back(chem_str_to_id.at(chem_str));
            reaction_ptr->reactant_deltas.emplace_back((unsigned int) std::stoul(coeff_st));
        }
    }

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
    if (r_filename.length() < 2) {
        std::cerr << "Error: Input file name " << r_filename << " is invalid." << std::endl;
        return 1;
    } else if (r_filename.substr(r_filename.length() - 2, 2) != ".r") {
        std::cerr << "Error: Input file type of " << r_filename << " is invalid." << std::endl;
        return 1;
    }

    std::ifstream r_file(r_filename);
    if (!r_file.is_open()) {
        std::cerr << "Error: Input .r file failed to be opened" << std::endl;
        return 1;
    }

    std::regex term_re("([A-Za-z'][A-Za-z0-9.'_]*) ([1-9][0-9]*)");
    std::regex field_re("^([^: ]+ [0-9]+[ ]*)+$");

    std::regex cr_char_re("\\r");

    unsigned int reaction_no = 0;
    std::vector<std::string> reaction_tokens;
    reaction_tokens.reserve(3);

    std::string temp("");
    while (getline(r_file, temp)) {
        temp = std::regex_replace(temp, cr_char_re, "");

        unsigned int num_splits = split_by(reaction_tokens, temp, ':');

        if (num_splits != 2) {
            print_error_msg_at_line(r_filename, reaction_no, "Missing field(s) in reaction.");
            return 1;
        }

        for (unsigned int i = 0; i < reaction_tokens.size(); ++i) {
            trim_left(reaction_tokens[i]);
            trim_right(reaction_tokens[i]);
        }

        if (reaction_tokens[RATE_FIELD].length() <= 0) {
            print_error_msg_at_line(r_filename, reaction_no, "Reaction rate is missing.");
            return 1;
        }

        try {
            float new_rate = std::stof(reaction_tokens[RATE_FIELD]);
            if (new_rate <= 0.0f) {
                print_error_msg_at_line(r_filename, reaction_no, "Reaction rate must be a positive non-zero value.");
                return 1;
            }

            crn.reactions.emplace_back(std::make_unique<reaction_t>());
            crn.reactions[reaction_no]->rate = new_rate;
        } catch (...) {
            print_error_msg_at_line(r_filename, reaction_no, "Reaction rate must be a positive non-zero value.");
            return 1;
        }

        std::smatch sm;
        if (!std::regex_match(reaction_tokens[REACTANT_FIELD], sm, field_re)) {
            print_error_msg_at_line(r_filename, reaction_no, "Term in reactant side of reaction is formatted incorrectly.");
            return 1;
        }

        if (fill_reaction_field(crn.reactions[reaction_no], chem_str_to_id, reaction_tokens[REACTANT_FIELD], REACTANT_FIELD, r_filename, reaction_no)) {
            return 1;
        } else if (fill_reaction_field(crn.reactions[reaction_no], chem_str_to_id, reaction_tokens[PRODUCT_FIELD], PRODUCT_FIELD, r_filename, reaction_no)) {
            return 1;
        }

        reaction_tokens.clear();
        reaction_no++;
    }

    r_file.close();
    return 0;
}