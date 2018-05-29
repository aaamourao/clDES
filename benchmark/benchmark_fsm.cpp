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
    QSet<cldes::ScalarType> non_contr;

    for (auto i = 0ul; i <= 100ul; ++i) {
        if (i % 2 == 0ul) {
            non_contr.insert(i);
        }
    }

    std::set<cldes::cldes_size_t> marked_states;
    marked_states.emplace(0);

    std::cout << "Generating Plants" << std::endl;

    cldes::DESystem c_1{2, 0, marked_states};
    c_1(0, 1) = 11ul;
    c_1(1, 0) = 12ul;

    cldes::DESystem c_2{2, 0, marked_states};
    c_2(0, 1) = 21ul;
    c_2(1, 0) = 22ul;

    cldes::DESystem milling{2, 0, marked_states};
    milling(0, 1) = 41ul;
    milling(1, 0) = 42ul;

    cldes::DESystem mp{2, 0, marked_states};
    mp(0, 1) = 81ul;
    mp(1, 0) = 82ul;

    cldes::DESystem lathe{3, 0, marked_states};
    lathe(0, 1) = 51ul;
    lathe(1, 0) = 52ul;
    lathe(0, 2) = 53ul;
    lathe(2, 0) = 54ul;

    cldes::DESystem c_3{3, 0, marked_states};
    c_3(0, 1) = 71ul;
    c_3(1, 0) = 72ul;
    c_3(0, 2) = 73ul;
    c_3(2, 0) = 74ul;

    cldes::DESystem robot{6, 0, marked_states};
    robot(0, 1) = 31ul;
    robot(1, 0) = 32ul;
    robot(0, 2) = 33ul;
    robot(2, 0) = 34ul;
    robot(0, 3) = 35ul;
    robot(3, 0) = 36ul;
    robot(0, 4) = 37ul;
    robot(4, 0) = 38ul;
    robot(0, 5) = 39ul;
    robot(5, 0) = 30ul;

    cldes::DESystem mm{4, 0, marked_states};
    mm(0, 1) = 61ul;
    mm(1, 2) = 63ul;
    mm(1, 3) = 65ul;
    mm(2, 0) = 64ul;
    mm(3, 0) = 66ul;

    std::cout << "Generating Specs" << std::endl;

    cldes::DESystem e_1{2, 0, marked_states};
    e_1(0, 1) = 12ul;
    e_1(1, 0) = 31ul;

    cldes::DESystem e_2{2, 0, marked_states};
    e_2(0, 1) = 22ul;
    e_2(1, 0) = 33ul;

    cldes::DESystem e_3{3, 0, marked_states};
    e_3(0, 1) = 32ul;
    e_3(1, 0) = 41ul;
    e_3(0, 2) = 42ul;
    e_3(2, 0) = 35ul;

    cldes::DESystem e_4{4, 0, marked_states};
    e_4(0, 1) = 34ul;
    e_4(1, 0) = 51ul;
    e_4(1, 0) = 53ul;
    e_4(0, 2) = 52ul;
    e_4(2, 0) = 37ul;
    e_4(0, 3) = 54ul;
    e_4(3, 0) = 39ul;

    cldes::DESystem e_5{2, 0, marked_states};
    e_5(0, 1) = 36ul;
    e_5(1, 0) = 61ul;

    cldes::DESystem e_6{2, 0, marked_states};
    e_6(0, 1) = 38ul;
    e_6(1, 0) = 63ul;

    cldes::DESystem e_7{3, 0, marked_states};
    e_7(0, 1) = 30ul;
    e_7(1, 0) = 71ul;
    e_7(0, 2) = 74ul;
    e_7(2, 0) = 65ul;

    cldes::DESystem e_8{3, 0, marked_states};
    e_8(0, 1) = 72ul;
    e_8(1, 0) = 81ul;
    e_8(0, 2) = 82ul;
    e_8(2, 0) = 73ul;

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
    std::cout << "Number of transitions of the plant " << plant.GetGraph().nnz()
              << std::endl;
    std::cout << "Computing the supervisor" << std::endl;
    std::cout << "Number of states of the spec: " << spec.Size() << std::endl;
    std::cout << "Number of transitions of the spec " << spec.GetGraph().nnz()
              << std::endl
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
              << supervisor.GetGraph().nnz() << std::endl;
}
