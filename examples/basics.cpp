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

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>
#include <set>
#include "cldes.hpp"

namespace ublas = boost::numeric::ublas;

int main() {
    const int n_states = 4;

    cldes::DESystem::StatesSet marked_states;
    marked_states.insert(0);
    marked_states.insert(2);

    const int init_state = 0;

    cldes::DESystem sys{n_states, init_state, marked_states};

    // Declare transitions: represented by prime numbers
    // TODO: Transitions and Events still need to be full implemented
    const float a = 2.0f;
    const float b = 3.0f;
    const float g = 5.0f;

    sys(0, 0) = a;
    sys(0, 2) = g;
    sys(1, 0) = a;
    sys(1, 1) = b;
    sys(2, 1) = a * g;
    sys(2, 2) = b;
    sys(2, 3) = a;
    sys(3, 1) = a;

    auto graph = sys.GetGraph();

    // std::cout << "Graph data: " << std::endl;
    // std::cout << graph << std::endl;

    if (&graph != NULL) {
        std::cout << "It is a beginning..." << std::endl;
    }

    auto accessible_states = sys.AccessiblePart();

    std::cout << "Accessible part: ";
    for (auto it = accessible_states.begin(); it != accessible_states.end();
         ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    return 0;
}
