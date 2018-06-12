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

 File: operations.hpp
 Description: Declaration of operation functions.
 =========================================================================
*/

#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "cldes/Constants.hpp"
#include <Eigen/Sparse>
#include <stack>
#include <tuple>

namespace cldes {

namespace op {
/*! \brief Graph represented by adjacency matrix type
 *
 * Graph of events represented by a bit array.
 */
template<uint8_t NEvents>
using GraphType = Eigen::SparseMatrix<EventsSet<NEvents>, Eigen::RowMajor>;

/*! \brief Vector of states type
 */
template<typename StorageIndex>
using StatesArray = std::vector<StorageIndex>;

/*! \brief tuple representing a state of a virtual synch (stage 1)
 *
 * (state_id_g0, state_id_g1)
 */
template<typename StorageIndex>
using StatesTupleHost = std::pair<StorageIndex, StorageIndex>;

/*! \brief Hash set of virtual states (stage 1)
 *
 * st = state_id_g1 * g0.size() + state_id_g0
 */
template<typename StorageIndex>
using StatesTableHost = spp::sparse_hash_set<StorageIndex>;

/*! \brief Hash map type for maps a virtual state to its new index
 *
 * It is necessary when states are removed.
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

/*! \brief Returns a system which represents a parallel composition
 * between two systems.
 *
 * The composed states are sorted by the right operand indexes:
 * e.g. indexes(sys0.size{3} || sys.size{2}) =
 * {0 = (0, 0), 1 = (1, 0), 2 = (2, 0), 3 = (0, 1), 4 = (1, 1), 5 = (2, 1)}
 *
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 */
template<uint8_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
Synchronize(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
            cldes::DESystem<NEvents, StorageIndex> const& aSys1);

/*! \brief Returns a virtual system which represents a parallel composition
 * between two systems.
 *
 * The composed states are sorted by the right operand indexes:
 * e.g. virtualsync(sys0.size{3} || sys.size{2}) =
 * {(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)}
 *
 * @param aVirtualSys Reference to the system which will be transformed.
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 */
template<uint8_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SynchronizeStage1(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

/*! \brief Transform a virtual system in a real system
 *
 * Calculate transitions and other parameters of a virtual system which
 * represents a virtual parallel composition.
 *
 * @param aVirtualSys Reference to the system which will be transformed.
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 */
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeStage2(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

/*! \brief Returns the output of the transitions function of a virtual system
 *
 * It assumes that the transition exists. It has a non defined behaviour when
 * the transition does not exist. It should be checked by the user before
 * calling TransitionVirtual.
 *
 * Implements:
 * f(q, e) = q_out;
 *
 * @param aP Plant system const reference
 * @param aE Spec system const reference
 * @param aQ The "from state"
 * @param aEvent Event of transition
 */
template<uint8_t NEvents, typename StorageIndex>
StorageIndex
TransitionVirtual(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1,
                  StorageIndex const& aQ,
                  cldes::ScalarType const& aEvent);

/*! \brief Remove bad states recursively
 *
 * Remove a state and all the states that become a bad state when the previous
 * one is removed.
 *
 * @param aVirtualSys reference to the virtual system which will have a state
 * removed
 * @param aP Plant system const reference
 * @param aE Specs system const reference
 * @param aNonContr Hash table containing all non-controllable events indexes.
 * @param aInvGraphP The plant bit array graph transposed.
 * @param aInvGraphE The spec bit array graph transposed.
 * @param aC A hash table containing the states currently added to the
 * virtual system.
 * @param aQ The state to be removed.
 * @param aNonContrBit Bit array containing the virtual system non controllable
 * events.
 * @param aRmTable A hash table containing all the removed states so far
 */
template<uint8_t NEvents, typename StorageIndex>
void
RemoveBadStates(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aE,
                GraphType<NEvents> const& aInvGraphP,
                GraphType<NEvents> const& aInvGraphE,
                StatesTableHost<StorageIndex>& aC,
                StorageIndex const& aQ,
                EventsSet<NEvents> const& aNonContrBit,
                StatesTableHost<StorageIndex>& aRmTable);

/*! \brief Returns the monolithic supervisor of a plant and a pec
 *
 * @param aP Plant system const reference
 * @param aE Specs system const reference
 * @param aNonContr Hash table containing all non-controllable events indexes.
 */
template<uint8_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SupervisorSynth(cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aE,
                EventsTableHost const& aNonContr);

} // namespace op
} // namespace cldes

// Including implementation
#include "cldes/src/operations/OperationsCore.hpp"
#endif // DESYSTEM_HPP
