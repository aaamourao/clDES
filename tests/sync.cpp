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

#include "cldes.hpp"
#include "operations/operations.hpp"

int main() {
    // Declare transitions: represented by prime numbers
    const float a = 2.0f;
    const float b = 3.0f;
    const float g = 5.0f;

    // Declare system G1
    int const nstatesG1 = 3;

    cldes::DESystem::StatesSet markedstatesG1;
    markedstatesG1.insert(0);
    markedstatesG1.insert(2);

    int const initstateG1 = 0;
    cldes::DESystem g1{nstatesG1, initstateG1, markedstatesG1};

    g1(0, 0) = a;
    g1(0, 2) = g;
    g1(1, 0) = a;
    g1(1, 1) = b;
    g1(2, 1) = a * g;
    g1(2, 2) = b;

    // Declare system G2
    int const nstatesG2 = 2;

    cldes::DESystem::StatesSet markedstatesG2;
    markedstatesG1.insert(1);

    int const initstateG2 = 0;
    cldes::DESystem g2{nstatesG2, initstateG2, markedstatesG2};

    g2(0, 0) = b;
    g2(0, 1) = a;
    g2(1, 0) = b;
    g2(1, 1) = a;

    auto stage1 = cldes::op::SynchronizeStage1(g1, g2);

    for (int i = 0; i < nstatesG1*nstatesG2; ++i) {
        std::cout << "(" << stage1[i].x0 << ", " << stage1[i].x1 << ")"
                  << std::endl;
    }

    std::cout << "Finishing test" << std::endl;

    return 0;
}
