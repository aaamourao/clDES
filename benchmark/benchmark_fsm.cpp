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
#include <cstdlib>
#include <vector>
#include "des/desystem.hpp"
#include "des/transition_proxy.hpp"
#include "operations/operations.hpp"
#include "testlib.hpp"

using namespace std::chrono;

int main() {
    std::set<cldes::cldes_size_t> marked_states;
    marked_states.emplace(0);

    QSet<cldes::ScalarType> non_contr;

    std::cout << "Generating Plants" << std::endl;

    cldes::DESystem c_1{2, 0, marked_states};
    c_1(0, 1) = 0ul;
    c_1(1, 0) = 1ul;

    non_contr.insert(1ul);

    cldes::DESystem c_2{2, 0, marked_states};
    c_2(0, 1) = 3ul;
    c_2(1, 0) = 4ul;

    non_contr.insert(4ul);

    cldes::DESystem milling{2, 0, marked_states};
    milling(0, 1) = 5ul;
    milling(1, 0) = 6ul;

    non_contr.insert(6ul);

    cldes::DESystem mp{2, 0, marked_states};
    mp(0, 1) = 7ul;
    mp(1, 0) = 8ul;

    non_contr.insert(8ul);

    cldes::DESystem lathe{3, 0, marked_states};
    lathe(0, 1) = 9ul;
    lathe(1, 0) = 10ul;
    lathe(0, 2) = 11ul;
    lathe(2, 0) = 12ul;

    non_contr.insert(10ul);
    non_contr.insert(12ul);

    cldes::DESystem c_3{3, 0, marked_states};
    c_3(0, 1) = 13ul;
    c_3(1, 0) = 14ul;
    c_3(0, 2) = 15ul;
    c_3(2, 0) = 16ul;

    non_contr.insert(14ul);
    non_contr.insert(16ul);

    cldes::DESystem robot{6, 0, marked_states};
    robot(0, 1) = 17ul;
    robot(1, 0) = 18ul;
    robot(0, 2) = 19ul;
    robot(2, 0) = 20ul;
    robot(0, 3) = 21ul;
    robot(3, 0) = 22ul;
    robot(0, 4) = 23ul;
    robot(4, 0) = 24ul;
    robot(0, 5) = 25ul;
    robot(5, 0) = 26ul;

    non_contr.insert(18ul);
    non_contr.insert(20ul);
    non_contr.insert(22ul);
    non_contr.insert(24ul);
    non_contr.insert(26ul);

    cldes::DESystem mm{4, 0, marked_states};
    mm(0, 1) = 27ul;
    mm(1, 2) = 28ul;
    mm(1, 3) = 29ul;
    mm(2, 0) = 30ul;
    mm(3, 0) = 31ul;

    non_contr.insert(30ul);
    non_contr.insert(31ul);

    std::cout << "Generating Specs" << std::endl;

    cldes::DESystem e_1{2, 0, marked_states};
    e_1(0, 1) = 1ul;
    e_1(1, 0) = 17ul;

    cldes::DESystem e_2{2, 0, marked_states};
    e_2(0, 1) = 4ul;
    e_2(1, 0) = 19ul;

    cldes::DESystem e_3{3, 0, marked_states};
    e_3(0, 1) = 18ul;
    e_3(1, 0) = 5ul;
    e_3(0, 2) = 6ul;
    e_3(2, 0) = 21ul;

    cldes::DESystem e_4{4, 0, marked_states};
    e_4(0, 1) = 20ul;
    e_4(1, 0) = 9ul;
    e_4(1, 0) = 11ul;
    e_4(0, 2) = 10ul;
    e_4(2, 0) = 23ul;
    e_4(0, 3) = 12ul;
    e_4(3, 0) = 25ul;

    cldes::DESystem e_5{2, 0, marked_states};
    e_5(0, 1) = 22ul;
    e_5(1, 0) = 27ul;

    cldes::DESystem e_6{2, 0, marked_states};
    e_6(0, 1) = 24ul;
    e_6(1, 0) = 28ul;

    cldes::DESystem e_7{3, 0, marked_states};
    e_7(0, 1) = 26ul;
    e_7(1, 0) = 13ul;
    e_7(0, 2) = 16ul;
    e_7(2, 0) = 29ul;

    cldes::DESystem e_8{3, 0, marked_states};
    e_8(0, 1) = 14ul;
    e_8(1, 0) = 7ul;
    e_8(0, 2) = 8ul;
    e_8(2, 0) = 15ul;

    std::vector<cldes::DESystem> plants = {c_1, c_2,   milling, lathe,
                                           mm,  robot, c_3,     mp};
    std::vector<cldes::DESystem> specs = {e_1, e_2, e_3, e_4,
                                          e_5, e_6, e_7, e_8};

    std::cout << "Synchronizing plants" << std::endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto last_result = plants[0ul];
    auto plant = last_result;
    PrintGraph(plants[0ul].GetGraph(), "plants[0]");
    for (auto i = 1ul; i < plants.size(); ++i) {
        plant = cldes::op::Synchronize(last_result, plants[i]);
        last_result = plant;
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Plants sync time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Synchronizing specs" << std::endl;
    t1 = high_resolution_clock::now();
    last_result = specs[0ul];
    auto spec = last_result;
    for (auto i = 1ul; i < specs.size(); ++i) {
        spec = cldes::op::Synchronize(last_result, specs[i]);
        last_result = spec;
    }
    t2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Specs sync time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << std::endl
              << "Number of states of plant: " << plant.Size() << std::endl;
    std::cout << "Number of transitions of the plant "
              << plant.GetGraph().nonZeros() << std::endl;
    std::cout << "Computing the supervisor" << std::endl;
    std::cout << "Number of states of the spec: " << spec.Size() << std::endl;
    std::cout << "Number of transitions of the spec "
              << spec.GetGraph().nonZeros() << std::endl
              << std::endl;

    t1 = high_resolution_clock::now();
    auto supervisor = cldes::op::SupervisorSynth(plant, spec, non_contr);
    t2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Supervisor synth time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Number of states of the supervisor: " << supervisor.Size()
              << std::endl;
    std::cout << "Number of transitions of the supervisor "
              << supervisor.GetGraph().nonZeros() << std::endl;
}
