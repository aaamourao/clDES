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

 File: examples/basics.cpp
 Description: Exemplify the basic usage of the library.
 =========================================================================
*/

#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "cldes.hpp"

namespace ublas = boost::numeric::ublas;

int main() {
    clDES::DESystem *sys;

    const int n_states = 6;

    const auto system_graph =
        new ublas::compressed_matrix<clDES::ScalarType>(n_states, n_states);

    // Declare transitions: represented by prime numbers
    // TODO: Transitions still need to be full implemented
    const float a = 2.0f;
    const float b = 3.0f;
    const float g = 5.0f;

    std::vector<int> marked_states;
    marked_states.push_back(0);
    marked_states.push_back(2);

    const int init_state = 0;

    (*system_graph)(0, 0) = a; (*system_graph)(0, 2) = g;
    (*system_graph)(1, 0) = a; (*system_graph)(1, 1) = b;
    (*system_graph)(2, 1) = a * g; (*system_graph)(2, 2) = b;

    sys = new clDES::DESystem(*system_graph, n_states, init_state, marked_states, false);

    ublas::compressed_matrix<clDES::ScalarType> graph = sys->GetGraph();

    if (&graph != NULL) {
        std::cout << "It is a beginning..." << std::endl;
    }

    delete sys;

    return 0;
}
