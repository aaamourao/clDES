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

 LacSED - Laboratorio de Analise e Controle de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: cldes/src/des/DESystemFwd.cpp
 Description: Forward declarations, includes and alias definitions.
 =========================================================================
*/
#ifndef DESYSTEM_CORE_HPP
#define DESYSTEM_CORE_HPP

#include <algorithm>
#include <Eigen/Sparse>
#include <boost/iterator/counting_iterator.hpp>
#include <functional>
#include <sparsepp/spp.h>
#include <stack>
#include <tuple>
#include <vector>

namespace cldes {

template<uint8_t NEvents, typename StorageIndex>
class DESystemCL;

/*! \brief Vector that contains arguments of inverse transition function
 *
 * V f(s_from, e) = s_to |-> f^-1(s_to, e) = s_from
 * f^-1 args are (s_to, e)
 */
template<typename StorageIndex>
using InvArgtrans = std::vector<std::pair<StorageIndex, cldes::ScalarType>>;

/*! \brief Hash map containing transitions
 *
 * V f(s_from, e) = s_to |-> f^-1(s_to, e) = s_from
 * {key = s_from, value= vec(s_from, e))
 */
template<typename StorageIndex>
using transMap = spp::sparse_hash_map<StorageIndex, InvArgtrans<StorageIndex>*>;

/*
 * Forward declarion of DESystemBase's friends class TransitionProxy. A
 * transition is an element of the adjascency matrix which implements
 * the des graph.
 */
template<uint8_t NEvents, typename StorageIndex>
class TransitionProxy;

/*
 * Forward declarion of DESystem class necessary for the forward declaration of
 * the DESystem's friend function op::synchronize
 */
template<uint8_t NEvents = kDefaultEventsN, typename StorageIndex = uint32_t>
class DESystem;

/*! \brief Vector of DES systems on host mem
 */
template<uint8_t NEvents, typename StorageIndex>
using DESVector = std::vector<DESystem<NEvents, StorageIndex>>;

// Forward declartions of friends functions which implement des operations
namespace op {
/*
 * Forward declaration of the synchronize virtual proxy
 */
template<uint8_t NEvents, typename StorageIndex>
class SyncSysProxy;
}

/*! \brief Alias for graph 3-tuple
 *
 * (s_from, s_to, transition_events)
 */
template<uint8_t NEvents>
using Triplet = Eigen::Triplet<EventsSet<NEvents>>;

/*! \brief Alias for bit graph 3-tuple
 *
 * (s_from, s_to, true)
 */
using BitTriplet = Eigen::Triplet<bool>;
}

#endif
