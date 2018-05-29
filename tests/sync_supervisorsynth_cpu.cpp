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

 File: test/kernels.hpp
 Description: Test cldes::op::Synchronize function, the parallel
 composition implementation.
 =========================================================================
*/

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include "cldes.hpp"

#include "operations/operations.hpp"
#include "testlib.hpp"

using namespace std::chrono;

int main() {
    // Declare transitions: represented by prime numbers
    cldes::ScalarType const a0 = 0; // event value: 1
    cldes::ScalarType const a1 = 1; // event value: 2
    cldes::ScalarType const b0 = 2; // event value: 4
    cldes::ScalarType const b1 = 3; // event value: 8

    QSet<cldes::ScalarType> non_contr = {b0, b1};

    std::set<cldes::cldes_size_t> plant_marked_states = {0};

    cldes::DESystem g1{2, 0, plant_marked_states};

    g1(0, 1) = a0;
    g1(1, 0) = b0;

    cldes::DESystem g2{2, 0, plant_marked_states};

    g2(0, 1) = a1;
    g2(1, 0) = b1;

    auto plant = cldes::op::Synchronize(g1, g2);

    PrintGraph(plant.GetGraph(), "Plant");

    std::set<cldes::cldes_size_t> spec_marked_states = {0, 1};

    cldes::DESystem spec{2, 0, spec_marked_states};

    spec(0, 1) = b0;
    spec(1, 0) = a1;

    PrintGraph(spec.GetGraph(), "Spec");

    auto t1 = high_resolution_clock::now();
    auto supervisor = cldes::op::SupervisorSynth(plant, spec, non_contr);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Supervisor synth time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Number of states of the supervisor: " << supervisor.Size()
              << std::endl;
    std::cout << "Number of transitions of the supervisor "
              << supervisor.GetGraph().nonZeros() << std::endl;

    std::ostringstream expected_result;

    expected_result << "0 1 0 0 0 0 " << std::endl;
    expected_result << "0 0 0 0 4 0 " << std::endl;
    expected_result << "8 0 0 1 0 0 " << std::endl;
    expected_result << "0 8 0 0 0 4 " << std::endl;
    expected_result << "0 0 2 0 0 0 " << std::endl;
    expected_result << "0 0 0 0 8 0 " << std::endl;
    expected_result << ">" << std::endl;
    ProcessResult(supervisor.GetGraph(), "< Sync graph",
                  expected_result.str().c_str());
    std::cout << "Synchronize time: " << duration << " microseconds"
              << std::endl;

    return 0;
}
