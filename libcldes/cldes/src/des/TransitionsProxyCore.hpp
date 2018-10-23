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

 File: cldes/src/des/TransitionProxyCore.hpp
 Description: TransitionProxy class implementation.
 =========================================================================
*/

template<uint8_t NEvents, typename StorageIndex>
cldes::TransitionProxy<NEvents, StorageIndex>::TransitionProxy(
  cldes::DESystem<NEvents, StorageIndex>* const aSysPtr,
  StorageIndex const& aLin,
  StorageIndex const& aCol)
  : sys_ptr_{ aSysPtr }
  , lin_{ aLin }
  , col_{ aCol }
{
}

template<uint8_t NEvents, typename StorageIndex>
cldes::TransitionProxy<NEvents, StorageIndex>&
cldes::TransitionProxy<NEvents, StorageIndex>::operator=(
  cldes::ScalarType aEventPos)
{
    // Add transition to the system
    sys_ptr_->events_[aEventPos] = true;

    // Create a unsigned long long representing the event
    EventsSet<NEvents> const event_ull{ 1ul << aEventPos };

    // Add transition to the state events hash table
    sys_ptr_->states_events_[lin_] |= event_ull;

    // Add transition to the state events inverted hash table
    sys_ptr_->inv_states_events_[col_] |= event_ull;

    // Add transition to graph
    EventsSet<NEvents> const last_value = sys_ptr_->graph_.coeff(lin_, col_);
    sys_ptr_->graph_.coeffRef(lin_, col_) = last_value | event_ull;
    sys_ptr_->graph_.makeCompressed();

    sys_ptr_->is_cache_outdated_ = true;

    return *this;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::TransitionProxy<NEvents, StorageIndex>::operator EventsSet<NEvents>()
  const
{
    return (sys_ptr_->graph_)(lin_, col_);
}
