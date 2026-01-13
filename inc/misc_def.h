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
 * This file contains misc. definitions used by cuAleae. Threshold types were based on threshold "codes" from aleae.h of Aleae's source code.
*/

#ifndef __MISC_DEF_H
#define __MISC_DEF_H

unsigned int typedef chem_id_t;

enum threshold_types {
    THRESH_LT = 0,
    THRESH_LE = 1,
    THRESH_GE = 2,
    THRESH_GT = 3,
    THRESH_N = 4
};

#endif