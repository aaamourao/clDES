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
#include <unordered_set>
#include <vector>
#include "des/desystem.hpp"
#include "operations/operations.hpp"
#include "testlib.hpp"

using namespace std::chrono;

void ClusterTool(unsigned int const &aNClusters,
                 std::vector<cldes::DESystem> &aPlants,
                 std::vector<cldes::DESystem> &aSpecs,
                 std::unordered_set<cldes::ScalarType> &non_contr);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Wrong usage");
    }

    std::vector<cldes::DESystem> plants;
    std::vector<cldes::DESystem> specs;
    std::unordered_set<cldes::ScalarType> non_contr;

    std::cout << "Generating ClusterTool(" << std::atoi(argv[1]) << ")"
              << std::endl;
    ClusterTool(std::atoi(argv[1]), plants, specs, non_contr);

    std::cout << "Synchronizing plants" << std::endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto last_result = plants[0];
    auto plant = last_result;
    PrintGraph(plants[0].GetGraph(), "plants[0]");
    for (auto i = 1; i < plants.size(); ++i) {
        plant = cldes::op::Synchronize(last_result, plants[i]);
        last_result = plant;
    PrintGraph(plants[i].GetGraph(), "plants[i]");
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Plants sync time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Synchronizing specs" << std::endl;
    t1 = high_resolution_clock::now();
    last_result = specs[0];
    PrintGraph(specs[0].GetGraph(), "spec[0]");
    auto spec = last_result;
    for (auto i = 1; i < specs.size(); ++i) {
        spec = cldes::op::Synchronize(last_result, specs[i]);
        last_result = spec;
    PrintGraph(specs[i].GetGraph(), "spec[i]");
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
