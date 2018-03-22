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

    for (auto it1 = graph.begin1(); it1 != graph.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            std::cout << "(" << it1.index1() << ", " << it2.index2()
                      << ") = " << *it2 << std::endl;
        }
    }

    auto accessible_states = sys.AccessiblePart();

    std::cout << "Accessible part: ";
    for (auto state : accessible_states) {
        std::cout << state << " ";
    }
    std::cout << std::endl;

    // Ublas way of initialize DESystem
    std::cout << "Initializing DESystem with an ublas compressed matrix"
              << std::endl;
    ublas::compressed_matrix<float> host_graph(n_states, n_states);

    // This graph has no transition from the 3rd state to th 4th one.
    host_graph(0, 0) = a;
    host_graph(0, 2) = g;
    host_graph(1, 0) = a;
    host_graph(1, 1) = b;
    host_graph(2, 1) = a * g;
    host_graph(2, 2) = b;
    host_graph(3, 1) = a;

    cldes::DESystem ublas_sys{host_graph, n_states, init_state, marked_states};

    auto ublas_graph = ublas_sys.GetGraph();

    for (auto it1 = ublas_graph.begin1(); it1 != ublas_graph.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            std::cout << "(" << it1.index1() << ", " << it2.index2()
                      << ") = " << *it2 << std::endl;
        }
    }

    auto ublas_accessible_states = ublas_sys.AccessiblePart();

    std::cout << "Accessible part: ";
    for (auto state : ublas_accessible_states) {
        std::cout << state << " ";
    }
    std::cout << std::endl;

    return 0;
}
