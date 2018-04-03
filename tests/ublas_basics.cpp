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

 File: tests/ublas_basics.cpp
 Description: Exemplify the basic usage of the library by initializing
 the discrete-event systems with ublas matrices.
 =========================================================================
*/

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include "cldes.hpp"

#include "testlib.hpp"

namespace ublas = boost::numeric::ublas;

int main() {
    std::cout << "Creating DES with ublas matrix" << std::endl;
    const int n_states = 4;

    cldes::DESystem::StatesSet marked_states;
    marked_states.insert(0);
    marked_states.insert(2);

    const int init_state = 0;

    ublas::compressed_matrix<float> adjmtr(n_states, n_states);

    // Declare transitions: represented by prime numbers
    const float a = 2.0f;
    const float b = 3.0f;
    const float g = 5.0f;

    adjmtr(0, 0) = a;
    adjmtr(0, 2) = g;
    adjmtr(1, 0) = a;
    adjmtr(1, 1) = b;
    adjmtr(2, 1) = a * g;
    adjmtr(2, 2) = b;
    adjmtr(2, 3) = a;

    cldes::DESystem sys{adjmtr, n_states, init_state, marked_states};
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
    ProcessResult(accessible_states, "Accessible part", "0 1 2 3");

    auto coaccessible_states = sys.CoaccessiblePart();
    ProcessResult(coaccessible_states, "Coaccessible part", "0 1 2");

    auto trimsys = sys.Trim();
    auto trimgraph = trimsys.GetGraph();
    PrintGraph(trimgraph, "Trim(Sys)");

    std::cout << "Creating new system with ublas matrix" << std::endl;
    ublas::compressed_matrix<float> host_graph(n_states, n_states);

    // This graph has no transition from the 3rd state to th 4th one.
    host_graph(0, 0) = a;
    host_graph(0, 2) = g;
    host_graph(1, 1) = b;
    host_graph(2, 1) = a * g;
    host_graph(2, 2) = b;
    host_graph(3, 1) = a;
    host_graph(3, 2) = a;

    cldes::DESystem ublas_sys{host_graph, n_states, init_state, marked_states};

    auto ublas_graph = ublas_sys.GetGraph();

    for (auto it1 = ublas_graph.begin1(); it1 != ublas_graph.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            std::cout << "(" << it1.index1() << ", " << it2.index2()
                      << ") = " << *it2 << std::endl;
        }
    }

    auto ublas_accessible_states = ublas_sys.AccessiblePart();
    ProcessResult(ublas_accessible_states, "Accessible part", "0 1 2");

    auto ublas_coaccessible_states = ublas_sys.CoaccessiblePart();
    ProcessResult(ublas_coaccessible_states, "Coaccessible part", "0 2 3");

    auto ublas_trimsys = ublas_sys.Trim();
    auto ublas_trimgraph = ublas_trimsys.GetGraph();
    PrintGraph(ublas_trimgraph, "Trim(New Sys)");

    return 0;
}
