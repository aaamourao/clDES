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

 File: desystem.cpp
 Description: DESystem class implementation. DESystem is a graph, which
 is modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#include "des/desystem.hpp"
#include <memory>
#include "viennacl/linalg/prod.hpp"

#include <iostream>

using namespace cldes;

DESystem::DESystem(GraphHostData &aGraph, int const &aStatesNumber,
                   int const &aInitState, std::vector<int> &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : graph_(new GraphHostData(aGraph)), init_state_(aInitState) {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        CacheGraph_();
    } else {
        device_graph_ = nullptr;
    }
}

DESystem::DESystem(int const &aStatesNumber, int const &aInitState,
                   std::vector<int> &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : graph_(new GraphHostData(aStatesNumber, aStatesNumber)),
      init_state_(aInitState) {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        CacheGraph_();
    } else {
        device_graph_ = nullptr;
    }
}

DESystem::~DESystem() {
    // Delete uBlas data
    if (graph_) {
        delete graph_;
    }
    if (dev_cache_enabled_) {
        delete device_graph_;
    }
}

DESystem::GraphHostData DESystem::GetGraph() const { return *graph_; }

std::set<int> DESystem::AccessiblePart() {
    // Cache graph temporally
    if (!dev_cache_enabled_) {
        CacheGraph_();
    } else if (is_cache_outdated_) {
        UpdateGraphCache_();
    }

    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */

    // Results sparse vector
    // TODO: When utils::Sum is implemented, change it to
    // viennacl::compressed_matrix
    ublas::vector<ScalarType> bfs_host_vector(states_number_);
    viennacl::vector<ScalarType> bfs_res_vector(states_number_);

    // Initialize sparse result vector
    for (int i = 0; i < states_number_; ++i) {
        bfs_host_vector[i] = 0;
    }
    bfs_host_vector(init_state_) = 1;

    // Copy initial vector to device memory
    viennacl::copy(bfs_host_vector, bfs_res_vector);

    viennacl::vector<ScalarType> X(states_number_);
    // Executes BFS
    for (int i = 0; i < states_number_; ++i) {
        if (i == 0) {
            X = bfs_res_vector;
        }
        auto Y = viennacl::linalg::prod(*device_graph_, X);
        bfs_res_vector += Y;
        X = Y;
    }

    // Remove graph_ from device memory, if it is set so
    if (!dev_cache_enabled_) {
        delete device_graph_;
    }

    // Copy result vector to host memory
    viennacl::copy(bfs_res_vector, bfs_host_vector);

    std::set<int> accessible_states;
    for (int i = 0; i < states_number_; ++i) {
        if (bfs_host_vector(i) != 0) {
            accessible_states.insert(i);
        }
    }

    return accessible_states;
}

void DESystem::CacheGraph_() {
    // If device graph is not allocated, allocate space for it
    if (!device_graph_) {
        device_graph_ = new viennacl::compressed_matrix<ScalarType>(
            states_number_, states_number_);
    }
    viennacl::copy(trans(*graph_), *device_graph_);

    is_cache_outdated_ = false;
}

// TODO: Is it always changing the value?
// So, in some situations, it may not be necessary to set
// is_cache_outdated_ = true
DESystem::GraphHostData::reference DESystem::operator()(int const &lin,
                                                        int const &col) {
    is_cache_outdated_ = true;

    return (*graph_)(lin, col);
}

void DESystem::UpdateGraphCache_() {
    viennacl::copy(trans(*graph_), *device_graph_);
    is_cache_outdated_ = false;
}
