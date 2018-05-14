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

 File: include/backend/kernels.cl
 Description: OpenCL kernels for the custom clDES operations.
 =========================================================================
*/

typedef struct StatesTuple {
    unsigned int x0;
    unsigned int x1;
} StatesTuple;

__kernel void Synchronize_Stage1(__global StatesTuple* aCTuples,
                                 unsigned int aBNumberStates) {
    // Workitems gets its index within index space
    int ix0 = get_global_id(0);
    int ix1 = get_global_id(1);

    unsigned int index = ix1 * aBNumberStates + ix0;

    aCTuples[index].x0 = ix0;
    aCTuples[index].x1 = ix1;
};
