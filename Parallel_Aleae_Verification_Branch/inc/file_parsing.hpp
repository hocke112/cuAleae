#ifndef __FILE_PARSING_HPP
#define __FILE_PARSING_HPP

#include <string>
#include <map>
#include "host_data_structures.hpp"

int parse_in_input_file(const std::string &in_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn);
int parse_r_input_file(const std::string &r_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn);

#endif