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
 * This file contains function prototypes for file parsing functions.
*/

#ifndef __FILE_PARSING_HPP
#define __FILE_PARSING_HPP

#include <string>
#include <map>
#include "host_data_structures.hpp"

int parse_in_input_file(const std::string &in_filename, std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn);
int parse_r_input_file(const std::string &r_filename, const std::map<std::string, chem_id_t> &chem_str_to_id, CRN &crn);

#endif