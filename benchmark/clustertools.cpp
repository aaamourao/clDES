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
        auto const istart = i * 8ul;

        if (i != aNClusters - 1ul) {
            cldes::DESystem r_i{4, 0, marked_states};
            r_i(0, 1) = istart;  // k
            r_i(1, 0) = istart + 1ul;  // k
            r_i(0, 2) = istart + 2ul;  // k
            r_i(2, 0) = istart + 3ul;  // k
            r_i(0, 3) = istart + 4ul;  // k
            r_i(3, 0) = istart + 5ul;  // k

            non_contr.insert(istart + 5ul);

            aPlants.push_back(r_i);
        } else {
            cldes::DESystem r_i{3, 0, marked_states};
            r_i(0, 1) = istart;  // k
            r_i(1, 0) = istart + 1ul;  // k
            r_i(0, 2) = istart + 4ul;  // k
            r_i(2, 0) = istart + 3ul;  // k

            aPlants.push_back(r_i);
        }

        non_contr.insert(istart + 1ul);
        non_contr.insert(istart + 3ul);

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


    for (auto i = 0ul; i < (aNClusters - 1); ++i) {
        auto const istart = 8ul * aNClusters + i * 4;

        cldes::DESystem e_ij{3, 0, marked_states};

        e_ij(0, 1) = istart;
        non_contr.insert(istart);

        e_ij(1, 0) = istart + 1ul;

        e_ij(0, 2) = istart + 2ul;
        non_contr.insert(istart + 2ul);

        e_ij(2, 0) = istart + 3ul;

        aSpecs.push_back(e_ij);
    }

    return;
}
