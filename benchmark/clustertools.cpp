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
#include "des/desystem.hpp"
#include "des/transition_proxy.hpp"

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
        for (auto k = 0ul; k <= 9ul; ++k) {
            auto index = i * 10ul + k;
            if (k % 2ul == 0ul) {
                non_contr.insert(index);
            }
        }
    }

    for (auto i = 0ul; i < aNClusters; ++i) {
        if (i != aNClusters - 1ul) {
            cldes::DESystem r_i{4, 0, marked_states};
            r_i(0, 1) = i * 10ul + 1ul;
            r_i(1, 0) = i * 10ul + 2ul;
            r_i(0, 2) = i * 10ul + 3ul;
            r_i(2, 0) = i * 10ul + 4ul;
            r_i(0, 3) = i * 10ul + 5ul;
            r_i(3, 0) = i * 10ul + 6ul;

            aPlants.push_back(r_i);
        } else {
            cldes::DESystem r_i{3, 0, marked_states};
            r_i(0, 1) = i * 10ul + 1ul;
            r_i(1, 0) = i * 10ul + 2ul;
            r_i(0, 2) = i * 10ul + 5ul;
            r_i(2, 0) = i * 10ul + 4ul;

            aPlants.push_back(r_i);
        }

        cldes::DESystem c_i{2, 0, marked_states};
        c_i(0, 1) = i * 10ul + 7ul;
        c_i(1, 0) = i * 10ul + 8ul;

        aPlants.push_back(c_i);

        cldes::DESystem e_i{3, 0, marked_states};
        e_i(0, 1) = i * 10ul + 2ul;
        e_i(1, 0) = i * 10ul + 7ul;
        e_i(0, 2) = i * 10ul + 8ul;
        e_i(2, 0) = i * 10ul + 5ul;

        aSpecs.push_back(e_i);
    }

    auto event_index_begin = 10 * aNClusters;

    for (auto i = 0ul; i < (aNClusters - 1); ++i) {
        cldes::DESystem e_ij{3, 0, marked_states};

        e_ij(0, 1) = event_index_begin + i * 4ul * aNClusters + 1ul;
        non_contr.insert(event_index_begin + i * 4ul * aNClusters + 1ul);

        e_ij(1, 0) = event_index_begin + i * 4ul * aNClusters + 2ul;

        e_ij(0, 2) = event_index_begin + i * 4ul * aNClusters + 3ul;
        non_contr.insert(event_index_begin + i * 4ul * aNClusters + 3ul);

        e_ij(2, 0) = event_index_begin + i * 4ul * aNClusters + 4ul;

        aSpecs.push_back(e_ij);
    }

    return;
}
