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

 File: examples/matrop.cpp
 Description: Exemplify how to use custom opencl matrix operations
 available on clDES.
 =========================================================================
*/

#include <iostream>
#include "backend/matrop.hpp"
#include "viennacl/vector.hpp"

int main() {
    int vector_size = 100;
    viennacl::vector<float> vec1(vector_size);
    viennacl::vector<float> vec2(vector_size);

    for (auto i = 0; i < vector_size; ++i) {
        vec1[i] = i;
        vec2[i] = i*2;
    }

    viennacl::result = cldes::backend::Sum(vec1, vec2);

    std::cout << "vec1 = " << vec1 << std::endl;
    std::cout << "vec2 = " << vec2 << std::endl;
    std::cout << "vec1 + vec2 = " << result << std::endl;

    return 0;
}
