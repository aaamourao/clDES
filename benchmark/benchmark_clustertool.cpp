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

void ClusterTool(unsigned long const &aNClusters,
                 std::vector<cldes::DESystem> &aPlants,
                 std::vector<cldes::DESystem> &aSpecs,
                 QSet<cldes::ScalarType> &non_contr) {
    if (aPlants.size() != 0 || aSpecs.size() != 0 || non_contr.size() ||
        aNClusters == 0) {
        throw std::runtime_error("ClusterTool: Invalid inputs");
    }

    std::set<cldes::cldes_size_t> marked_states;
    marked_states.emplace(0);

    for (auto i = 0ul; i < aNClusters; ++i) {
        auto const istart = i * 8ul;

        if (i != aNClusters - 1ul) {
            cldes::DESystem r_i{4, 0, marked_states};
            r_i(0, 1) = istart;        // k
            r_i(1, 0) = istart + 1ul;  // k
            r_i(0, 2) = istart + 2ul;  // k
            r_i(2, 0) = istart + 3ul;  // k
            r_i(0, 3) = istart + 4ul;  // k
            r_i(3, 0) = istart + 5ul;  // k

            aPlants.push_back(r_i);
        } else {
            cldes::DESystem r_i{3, 0, marked_states};
            r_i(0, 1) = istart;        // k
            r_i(1, 0) = istart + 1ul;  // k
            r_i(0, 2) = istart + 4ul;  // k
            r_i(2, 0) = istart + 3ul;  // k

            aPlants.push_back(r_i);
        }

        non_contr.insert(istart + 1ul);
        non_contr.insert(istart + 3ul);
        non_contr.insert(istart + 5ul);

        cldes::DESystem c_i{2, 0, marked_states};
        c_i(0, 1) = istart + 6ul;  // k
        c_i(1, 0) = istart + 7ul;  // k

        non_contr.insert(istart + 7ul);

        aPlants.push_back(c_i);

        cldes::DESystem e_i{3, 0, marked_states};
        e_i(0, 1) = istart + 1ul;  // k
        e_i(1, 0) = istart + 6ul;  // k
        e_i(0, 2) = istart + 7ul;  // k
        e_i(2, 0) = istart + 4ul;  // k

        aSpecs.push_back(e_i);
    }

    for (auto i = 1ul; i < aNClusters; ++i) {
        auto const istart = (i - 1) * 8ul;

        cldes::DESystem e_ij{3, 0, marked_states};

        e_ij(0, 1) = istart + 5ul;
        e_ij(1, 0) = istart + 8ul;
        e_ij(0, 2) = istart + 8ul + 3ul;
        e_ij(2, 0) = istart + 2ul;

        non_contr.insert(istart + 5ul);
        non_contr.insert(istart + 8ul + 3ul);

        aSpecs.push_back(e_ij);
    }

    return;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Wrong usage");
    }

    std::set<cldes::cldes_size_t> marked_states;
    cldes::DESystem plant{1, 0, marked_states};
    cldes::DESystem spec{1, 0, marked_states};
    QSet<cldes::ScalarType> non_contr;

    high_resolution_clock::time_point t1;
    high_resolution_clock::time_point t2;

    {
        std::vector<cldes::DESystem> plants;
        std::vector<cldes::DESystem> specs;

        std::cout << "Generating ClusterTool(" << std::atoi(argv[1]) << ")"
                  << std::endl;
        ClusterTool(std::atoi(argv[1]), plants, specs, non_contr);

        std::cout << "Synchronizing plants" << std::endl;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        auto last_result = plants[0ul];
        plant = last_result;
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
        spec = last_result;
        for (auto i = 1ul; i < specs.size(); ++i) {
            spec = cldes::op::Synchronize(last_result, specs[i]);
            last_result = spec;
        }
        t2 = high_resolution_clock::now();
        duration = duration_cast<microseconds>(t2 - t1).count();

        std::cout << "Specs sync time spent: " << duration << " microseconds"
                  << std::endl;
    }

    std::cout << std::endl
              << "Number of states of plant: " << plant.Size() << std::endl;
    std::cout << "Number of transitions of the plant "
              << plant.GetGraph().nonZeros() << std::endl;
    std::cout << "Number of states of the spec: " << spec.Size() << std::endl;
    std::cout << "Number of transitions of the spec "
              << spec.GetGraph().nonZeros() << std::endl
              << std::endl;

    std::cout << "{plant, spec}.Trim()" << std::endl;
    t1 = high_resolution_clock::now();
    plant.Trim();
    spec.Trim();
    t2 = high_resolution_clock::now();

    std::cout << std::endl
              << "Number of states of plant: " << plant.Size() << std::endl;
    std::cout << "Number of transitions of the plant "
              << plant.GetGraph().nonZeros() << std::endl;
    std::cout << "Number of states of the spec: " << spec.Size() << std::endl;
    std::cout << "Number of transitions of the spec "
              << spec.GetGraph().nonZeros() << std::endl
              << std::endl;

    auto duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Trim time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Computing the supervisor" << std::endl;
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
