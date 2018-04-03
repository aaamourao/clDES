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

 File: tests/basics.cpp
 Description: Exemplify the basic usage of the library.
 =========================================================================
*/

#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include "cldes.hpp"

template <typename T, typename StringType>
std::ostringstream ReadResult(T const &aOpResult, StringType const aHeader) {
    std::ostringstream result;

    result << aHeader << ": ";
    for (auto state : aOpResult) {
        result << state << " ";
    }
    result << std::endl;

    return result;
}

template <typename StringType>
void TestResult(StringType const &aResult, StringType const aExpected) {
    std::cout << aResult;
    assert(aResult.find(aExpected) != std::string::npos);
}

template <typename T, typename StringType>
void ProcessResult(T &aOpResult, StringType const aHeader,
                   StringType const aExpected) {
    auto result = ReadResult(aOpResult, aHeader);
    std::ostringstream expected;
    expected << aHeader << ": " << aExpected;
    TestResult(result.str(), expected.str());
}

int main() {
    std::cout << "Creating DES" << std::endl;
    const int n_states = 4;

    cldes::DESystem::StatesSet marked_states;
    marked_states.insert(0);
    marked_states.insert(2);

    const int init_state = 0;

    cldes::DESystem sys{n_states, init_state, marked_states};

    // Declare transitions: represented by prime numbers
    const float a = 2.0f;
    const float b = 3.0f;
    const float g = 5.0f;

    sys(0, 0) = a;
    sys(0, 2) = g;
    sys(1, 0) = a;
    sys(1, 1) = b;
    sys(2, 1) = a * g;
    sys(2, 2) = b;
    sys(2, 3) = a;

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

    sys.Trim();

    std::cout << "Creating new system" << std::endl;

    cldes::DESystem new_sys{n_states, init_state, marked_states};

    // This graph has no transition from the 3rd state to th 4th one.
    new_sys(0, 0) = a;
    new_sys(0, 2) = g;
    new_sys(1, 1) = b;
    new_sys(2, 1) = a * g;
    new_sys(2, 2) = b;
    new_sys(3, 1) = a;
    new_sys(3, 2) = a;

    auto new_graph = new_sys.GetGraph();

    for (auto it1 = new_graph.begin1(); it1 != new_graph.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            std::cout << "(" << it1.index1() << ", " << it2.index2()
                      << ") = " << *it2 << std::endl;
        }
    }

    auto new_accessible_states = new_sys.AccessiblePart();

    ProcessResult(new_accessible_states, "Accessible part", "0 1 2");

    auto new_coaccessible_states = new_sys.CoaccessiblePart();

    ProcessResult(new_coaccessible_states, "Coaccessible part", "0 2 3");

    new_sys.Trim();

    return 0;
}