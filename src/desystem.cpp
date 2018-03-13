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

DESystem::DESystem(ublas::compressed_matrix<ScalarType> &aGraph,
                   const int &aStatesNumber, const int &aInitState,
                   std::vector<int> &aMarkedStates,
                   const bool &aDevCacheEnabled)
    : graph_(&aGraph), init_state_(aInitState) {
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
    // TODO: Is device_graph_ on heap?
    if (dev_cache_enabled_) {
        delete device_graph_;
    }
}

ublas::compressed_matrix<ScalarType> DESystem::GetGraph() const {
    return *graph_;
}

std::set<int> DESystem::AccessiblePart() {
    // Cache graph temporally
    if (!dev_cache_enabled_) {
        CacheGraph_();
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
}

ublas::compressed_matrix<ScalarType>::reference DESystem::operator()(int &lin,
                                                                     int &col) {
    return (*graph_)(lin, col);
}
