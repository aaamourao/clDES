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

 File: cldes/src/operations/OperationsCore.cpp
 =========================================================================
*/

#ifndef OPERATIONS_FWD_HPP
#define OPERATIONS_FWD_HPP

#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <stack>
#include <tuple>

namespace cldes {
namespace op {
/*! \brief tuple representing a state of a virtual synch (stage 1)
 * \details (state_id_g0, state_id_g1)
 */
template<typename StorageIndex>
using StatesTupleHost = std::pair<StorageIndex, StorageIndex>;

/*! \brief Hash set of virtual states (stage 1)
 * \details st = state_id_g1 * g0.size() + state_id_g0
 */
template<typename StorageIndex>
using StatesTableHost = spp::sparse_hash_set<StorageIndex>;

/*! \brief Hash map type for maps a virtual state to its new index
 * \details It is necessary when states are removed.
 */
template<typename StorageIndex>
using SparseStatesMap = spp::sparse_hash_map<StorageIndex, StorageIndex>;

/*! \brief Stack of states type
 */
template<typename StorageIndex>
using StatesStack = std::stack<StorageIndex>;

/*! \brief Hash set containing events indexes
 */
using EventsTableHost = spp::sparse_hash_set<uint8_t>;

template<uint8_t NEvents, typename StorageIndex>
using Node = std::shared_ptr<DESystemBase<NEvents, StorageIndex, SyncSysProxy<NEvents, StorageIndex>>>;

template<uint8_t NEvents, typename StorageIndex>
using RefNodes = std::vector<Node<NEvents, StorageIndex>>;

template<uint8_t NEvents, typename StorageIndex>
using BinExprTree =
  std::pair<Node<NEvents, StorageIndex>, RefNodes<NEvents, StorageIndex>>;
}
}

#endif // OPERATIONS_FWD_HPP
