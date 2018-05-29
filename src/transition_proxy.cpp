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
        std::string error = "ValueError: Max transition value = " +
                            std::to_string(g_max_events);
        throw error;
    }

    // Add transition to the system
    sys_ptr_->events_[aEventPos] = true;

    // Create a unsigned long long representing the event
    DESystem::EventsSet event_ull;
    event_ull[aEventPos] = 1ull;

    // Add transition to the state events hash table
    if (sys_ptr_->states_events_.find(lin_) != sys_ptr_->states_events_.end()) {
        sys_ptr_->states_events_[lin_] =
            sys_ptr_->states_events_[lin_] | event_ull;
    } else {
        sys_ptr_->states_events_[lin_] = event_ull;
    }

    // Add transition to the state events inverted hash table
    if (sys_ptr_->inv_states_events_.find(col_) !=
        sys_ptr_->inv_states_events_.end()) {
        sys_ptr_->inv_states_events_[col_] =
            sys_ptr_->inv_states_events_[col_] | event_ull;
    } else {
        sys_ptr_->inv_states_events_[col_] = event_ull;
    }

    // Add transition to graph
    if ((sys_ptr_->graph_)(lin_, col_) == 0) {
        (sys_ptr_->graph_)(lin_, col_) = event_ull;
    } else {
        DESystem::EventsSet last_value = (sys_ptr_->graph_)(lin_, col_);
        (sys_ptr_->graph_)(lin_, col_) = last_value | event_ull;
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
