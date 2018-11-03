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

 LacSED - Laboratorio de Analise e Controle de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: test/kernels.hpp
 Description: Test cldes::op::Synchronize function, the parallel
 composition implementation.
 =========================================================================
*/

#define SPP_ALLOC_SZ 1

#include "cldes/DESystem.hpp"
#include "cldes/operations/Operations.hpp"
#include "cldes/operations/SuperProxy.hpp"
#include "clustertool.hpp"
#include "testlib.hpp"
#include <chrono>
#include <cstdlib>
#include <vector>

using namespace std::chrono;

int
main()
{
    using StorageIndex = uint64_t;

    std::set<StorageIndex> marked_states;
    cldes::DESystem<56, StorageIndex> plant{ 1, 0, marked_states };
    cldes::DESystem<56, StorageIndex> spec{ 1, 0, marked_states };
    cldes::DESystem<56>::EventsTable non_contr;

    high_resolution_clock::time_point t1;
    high_resolution_clock::time_point t2;

    {
        std::vector<cldes::DESystem<56, StorageIndex>> plants;
        std::vector<cldes::DESystem<56, StorageIndex>> specs;

        std::cout << "Generating ClusterTool(7)" << std::endl;
        ClusterTool(7, plants, specs, non_contr);
        std::cout << std::endl;

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
    auto supervisor =
      cldes::op::SuperProxy<56, StorageIndex>(plant, spec, non_contr);
    t2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(t2 - t1).count();
    std::cout << "Supervisor synth time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << std::endl;
    std::cout << "Number of states of the supervisor: " << supervisor.Size()
              << std::endl;
    // std::cout << "Number of transitions of the supervisor "
    //           << supervisor.GetGraph().nonZeros() << std::endl;

    return 0;
}