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

 File: cldes/src/operations/SyncSysProxyCore.hpp
 Description: includes and alias.
 =========================================================================
*/

#include <algorithm>

namespace cldes {
namespace op {

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeEmptyStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

// Alias to events hash map
using EventsTableHost = spp::sparse_hash_set<uint8_t>;

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(
  DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>> const&
    aP,
  DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>> const&
    aE,
  EventsTableHost const& aNonContr);
}
}
