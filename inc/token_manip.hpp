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
 * This file contains the logic to create and manipulate tokens.
*/

#ifndef __TOKEN_MANIP_HPP
#define __TOKEN_MANIP_HPP

#include <vector>
#include <string>

void trim_left(std::string &s);
void trim_right(std::string &s);
void trim(std::string &s);
void trim_tokens(std::vector<std::string> &tokens);

void split_by(std::vector<std::string> &tokens, const std::string &s, char delim);

#endif