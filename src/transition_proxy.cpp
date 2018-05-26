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

 File: transition_proxy.cpp
 Description: Proxy to an element of the graph_ data member from DESystem.
 =========================================================================
*/

#include "des/transition_proxy.hpp"
#include "des/desystem.hpp"

using namespace cldes;

TransitionProxy::TransitionProxy(DESystem *const aSysPtr,
                                 cldes_size_t const &aLin,
                                 cldes_size_t const &aCol)
    : sys_ptr_{aSysPtr}, lin_{aLin}, col_{aCol} {}

TransitionProxy &TransitionProxy::operator=(ScalarType aEventPos) {
    if (aEventPos > g_max_events) {
        throw "ValueError: Max transition value = 64";
    }

    // Add transition to the system and to the state
    sys_ptr_->events_[aEventPos] = true;
    (sys_ptr_->states_events_[lin_])[aEventPos] = true;

    // Add transition to graph
    std::bitset<64> event_l;
    event_l[aEventPos] = true;
    if ((sys_ptr_->graph_)(lin_, col_) == 0) {
        (sys_ptr_->graph_)(lin_, col_) = event_l;
    } else {
        std::bitset<64> last_value = ((sys_ptr_->graph_)(lin_, col_));
        (sys_ptr_->graph_)(lin_, col_) = last_value | event_l;
    }

    // Add transition to bit graph, which is transposed
    (sys_ptr_->bit_graph_)(col_, lin_) = true;

    sys_ptr_->is_cache_outdated_ = true;

    return *this;
}

/*
TransitionProxy::operator ScalarType() {
    return (sys_ptr_->graph_)(lin_, col_);
}
*/
