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
#include <set>
#include <vector>
#include "clustertools.hpp"
#include "des/desystem.hpp"

using namespace std::chrono;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Wrong usage");
    }

    std::vector<cldes::DESystem> plants;
    std::vector<cldes::DESystem> specs;
    std::set<cldes::ScalarType> non_contr;

    std::cout << "Generating ClusterTool(" << std::atoi(argv[1]) << ")"
              << std::endl;
    ClusterTool(std::atoi(argv[1]), plants, specs, non_contr);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    DESystem plant = *(plants.back());
    plants.pop_back();
    for (auto sys : plants) {
        DESystem op0 = plant;
        DESystem op1 = sys;
        plant = cldes::op::Synchronize(op0, op1);
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Plants sync time spent: " << duration << " microseconds"
              << std::endl;

    t1 = high_resolution_clock::now();
    DESystem spec = *(specs.back());
    specs.pop_back();
    for (auto sys : specs) {
        DESystem op0 = spec;
        DESystem op1 = sys;
        spec = cldes::op::Synchronize(op0, op1);
    }
    t2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Specs sync time spent: " << duration << " microseconds"
              << std::endl;

    t1 = high_resolution_clock::now();
    auto supervisor = cldes::op::SupervisorSynth(plants, specs, non_contr);
    t2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Supervisor synth time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Number of states of the supervisor: " << supervisor.Size()
              << std::endl;
    std::cout << "Number of transitions of the supervisor "
              << supervisor.GetGraph().nnz() << std::endl;
}
