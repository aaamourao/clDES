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

#include <vector>
#include "des/desystem.cpp"

std::vector<unsigned int> GeneratePrimeNumbers(unsigned int aMax) {
    std::vector<cldes::ScalarType> primes;
    primes.push_back(2);
    for (int i = 3; i < aMax; i++) {
        bool prime = true;
        for (auto j = 0; j < primes.size() && primes[j] * primes[j] <= i; j++) {
            if (i % primes[j] == 0) {
                prime = false;
                break;
            }
        }
        if (prime) {
            primes.push_back(static_cast<cldes::ScalarType>(i));
            cout << i << " ";
        }
    }
    return primes;
}

void ClusterTool(unsigned int const &aNClusters,
                 std::vector<cldes::DESystem> &aPlants,
                 std::vector<cldes::DESystem> &aSpecs,
                 std::set<cldes::DESystem> &non_contr) {
    if (aPlants.size() != 0 || aSpecs.size() != 0 || non_contr.size() ||
        aNClusters == 0) {
        throw std::runtime_error("ClusterTool: Invalid inputs");
    }

    std::set<cldes_size_t> marked_states;
    marked_states.emplace(0);

    auto n_events = 12 * aNClusters;
    auto events = GeneratePrimeNumbers(n_events);

    std::set<cldes::ScalarType> non_contr;
    for (auto i = 0; i < aNClusters; ++i) {
        for (auto k = 0; k < 9; ++k) {
            auto index = i * 9 + k;
            if (index < n_events && k % 2 != 0) {
                non_contr.emplace(events[index]);
            }
        }
    }

    for (auto i = 0; i < aNClusters; ++i) {
        if (i != aNClusters - 1) {
            cldes::DESystem r_i{4, 0, marked_states};
            r_i(0, 1) = events[i * aNClusters];
            r_i(1, 0) = events[i * aNClusters + 1];
            r_i(0, 2) = events[i * aNClusters + 2];
            r_i(2, 0) = events[i * aNClusters + 3];
            r_i(0, 3) = events[i * aNClusters + 4];
            r_i(3, 1) = events[i * aNClusters + 5];
            r_i.InsertEvents();
        } else {
            cldes::DESystem r_i{3, 0, marked_states};
            r_i(0, 1) = events[i * aNClusters];
            r_i(1, 0) = events[i * aNClusters + 1];
            r_i(0, 2) = events[i * aNClusters + 4];
            r_i(2, 0) = events[i * aNClusters + 3];
            r_i.InsertEvents();
        }

        aPlants.push_back(r_i);

        cldes::DESystem c_i{2, 0, marked_states};
        c_i(0, 1) = events[i * aNClusters + 6];
        c_i(1, 0) = events[i * aNClusters + 7];
        c_i.InsertEvents();

        aPlants.push_back(c_i);

        cldes::DESystem e_i{3, 0, marked_states};
        e_i(0, 1) = events[i * aNClusters + 1];
        e_i(1, 0) = events[i * aNClusters + 6];
        e_i(0, 2) = events[i * aNClusters + 7];
        e_i(2, 0) = events[i * aNClusters + 4];
        e_i.InsertEvents();

        aSpecs.push_back(e_i)
    }

    auto event_index_begin = 8 * aNClusters;

    for (auto i = 0; i < (aNClusters - 1); ++i) {
        cldes::DESystem e_ij{3, 0, marked_states};
        e_ij(0, 1) = events[event_index_begin + i * aNClusters + 1];
        e_ij(1, 0) = events[event_index_begin + i * aNClusters];
        e_ij(0, 2) = events[event_index_begin + i * aNClusters + 3];
        e_ij(2, 0) = events[event_index_begin + i * aNClusters + 2];
        e_ij.InsertEvents();

        aSpecs.push_back(e_i)
    }

    return;
}
