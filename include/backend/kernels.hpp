/*
 =========================================================================
 This file is part of clDES

 clDES: an OpenCL library for Discrete Event Systems computing.

 clDES is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 clDES is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with clDES.  If not, see <http://www.gnu.org/licenses/>.

 Copyright (c) 2018 - Adriano Mourao <adrianomourao@protonmail.com>
                      madc0ww @ [https://github.com/madc0ww]

 LacSED - Laboratório de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: include/backend/kernels.hpp
 Description: OpenCL kernels for the custom clDES operations.
 =========================================================================
*/

namespace cldes::backend {

char const* const cldes_kernels =
    "__kernel void elementwise_add(\n"
    "          __global const unsigned int * A_row_indices,\n"
    "          __global const unsigned int * A_col_indices, \n"
    "          __global const unsigned float * A_elements, \n"
    "          unsigned int * A_size, \n"
    "          __global const unsigned int * B_row_indices,\n"
    "          __global const unsigned int * B_col_indices, \n"
    "          __global const unsigned float * B_elements, \n"
    "          unsigned int * B_size) \n"
    "{ \n"
    // Workitems gets its index within index space
    "  int const ix = get_global_id(0); \n"
    "  int const iy = get_global_id(1); \n"
    "  for (unsigned int i = get_global_id(0); i < size; i += "
    "      get_global_size(0))\n"
    "    result[i] = vec1[i] + vec2[i];\n"
    "};\n";

}
