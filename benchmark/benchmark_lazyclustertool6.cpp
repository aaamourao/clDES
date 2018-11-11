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
 Description: Test cldes::op::synchronize function, the parallel
 composition implementation.
 =========================================================================
*/

#include "cldes/DESystem.hpp"
#include "cldes/operations/Operations.hpp"
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

    cldes::DESystem<48, StorageIndex>::EventsTable non_contr;
    cldes::DESVector<48, StorageIndex> plants;
    cldes::DESVector<48, StorageIndex> specs;

    std::cout << "Generating ClusterTool(6)" << std::endl;
    ClusterTool(6, plants, specs, non_contr);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    auto supervisor = cldes::op::supC(plants, specs, non_contr);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(t2 - t1).count();

    std::cout << std::endl;
    std::cout << "Supervisor synth time spent: " << duration << " microseconds"
              << std::endl;

    std::cout << "Number of states of the supervisor: " << supervisor.Size()
              << std::endl;
    std::cout << "Number of transitions of the supervisor "
              << supervisor.getGraph().nonZeros() << std::endl;
}
