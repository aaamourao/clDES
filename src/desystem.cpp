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

#include "cldes.hpp"
#include "viennacl/linalg/prod.hpp"

using namespace clDES;

DESystem::DESystem(ublas::compressed_matrix<ScalarType> *aGraph,
                   const int &aStatesNumber, const int &aInitState,
                   std::vector<int> aMarkedStates, const bool &aDevCacheEnabled)
    : graph_(aGraph), init_state_(aInitState) {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
}

DESystem::~DESystem() {
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

ublas::compressed_matrix<ScalarType> &DESystem::AccessiblePart() const {
    auto result_vector = new ublas::compressed_matrix<ScalarType>();

    return *result_vector;
}
