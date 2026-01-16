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
 * This file contains the function prototypes for creating and manipulating tokens.
*/

#include <cctype>

#include "token_manip.hpp"

/*
 * A helper function to remove whitespace at beginning of a string.
 * @param s: a string to be trimmed
*/
void trim_left(std::string &s) {
    if (s.length() < 1)
        return;

    unsigned int i = 0;
    while (i <= s.length() && std::isspace(s[i])) {
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
    while (i > 0 && std::isspace(s[i])) {
        --i;
    }

    if (s[i] != ' ')
        ++i;

    s.erase(i, s.length());
}

void trim(std::string &s) {
    trim_left(s);
    trim_right(s);
}

void trim_tokens(std::vector<std::string> &tokens) {
    for (unsigned int i = 0; i < tokens.size(); ++i) {
        trim(tokens[i]);
    }
}

/*
 * A helper function to split a string into a list of tokens.
 *
 * @param tokens: an empty list of tokens that will be filled
 * @param s: an input string to get tokens from
 * @char delim: a delimiting character to split the string
*/
void split_by(std::vector<std::string> &tokens, const std::string &s, char delim) {
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
}