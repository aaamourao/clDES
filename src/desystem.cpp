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

// TODO: delete this before commit
#include <iostream>

using namespace cldes;

DESystem::DESystem(GraphHostData const &aGraph,
                   cldes_size_t const &aStatesNumber,
                   cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : graph_{new GraphHostData{aGraph}}, init_state_{aInitState} {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;
    device_graph_ = nullptr;

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        CacheGraph_();
    }
}

DESystem::DESystem(cldes_size_t const &aStatesNumber,
                   cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : DESystem::DESystem{GraphHostData{aStatesNumber, aStatesNumber},
                         aStatesNumber, aInitState, aMarkedStates,
                         aDevCacheEnabled} {}

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

DESystem::StatesSet DESystem::AccessiblePart() {
    // Cache graph temporally
    if (!dev_cache_enabled_) {
        CacheGraph_();
    } else if (is_cache_outdated_) {
        UpdateGraphCache_();
    }

    // Executes a BFS on graph_
    auto accessible_states = Bfs_();

    // Remove graph_ from device memory, if it is set so
    if (!dev_cache_enabled_) {
        delete device_graph_;
    }

    // Avoiding memory leak
    // Workaround due to set issue on Bfs_
    auto ret_acc_states = StatesSet{*accessible_states};
    delete accessible_states;

    return ret_acc_states;
}

void DESystem::CacheGraph_() {
    // Allocate space for device_graph_
    if (device_graph_ == nullptr) {
        device_graph_ = new viennacl::compressed_matrix<ScalarType>{
            states_number_, states_number_};
    }
    viennacl::copy(trans(*graph_), *device_graph_);

    is_cache_outdated_ = false;
}

// TODO: Is it always changing the value?
// So, in some situations, it may not be necessary to set
// is_cache_outdated_ = true
DESystem::GraphHostData::reference DESystem::operator()(
    cldes_size_t const &lin, cldes_size_t const &col) {
    is_cache_outdated_ = true;

    return (*graph_)(lin, col);
}

void DESystem::UpdateGraphCache_() {
    viennacl::copy(trans(*graph_), *device_graph_);
    is_cache_outdated_ = false;
}

DESystem::StatesSet *DESystem::Bfs_(cldes_size_t const &aInitialNode) {
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    StatesVector host_x{states_number_, 1};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x(aInitialNode, 0) = 1;

    // Copy searching node to device memory
    StatesDeviceVector x{states_number_, 1};
    viennacl::copy(host_x, x);

    /*
     * Only Odin knows why I can't return by value here. When I do that
     * it seems that the compiler was returning by reference, but, of course,
     * the data were deleted, since it valid only on this scope.
     *
     * on this scope, gdb was poiting a 4 elements set
     * when calling the function a corupted set
     * $1 = std::set with 7138656 elements<error reading variable: Cannot access
     * memory at address 0xc99000000d44e8d7>
     *
     * segfault
     * 0x00007ffff7581903 in std::local_Rb_tree_increment
     * (__x=0xc99000000d44e8c7) at
     * ../../../../../libstdc++-v3/src/c++98/tree.cc:65
     * 65              while (__x->_M_left != 0) )
     */
    auto accessed_states = new StatesSet;

    // Executes BFS
    for (auto i = 0; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        StatesDeviceVector y = viennacl::linalg::prod(*device_graph_, x);
        x = y;

        // Unfortunatelly, until now, ViennaCL does not allow iterating on
        // compressed matrices. Until it is implemented, it is necessary
        // to copy the vector to the host memory.
        viennacl::copy(x, host_x);
        for (auto state = host_x.begin1(); state != host_x.end1(); ++state) {
            accessed_states->emplace(state.index1());
        }
    }

    return accessed_states;
}
