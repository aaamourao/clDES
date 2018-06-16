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

 File: cldes/operations/Operations.hpp
 Description: Definition of operation functions: Parallel composition,
 virtual parallel composition and supervisor synthesis.
 =========================================================================
*/

#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "cldes/Constants.hpp"
#include "cldes/DESystem.hpp"
#include "cldes/operations/SyncSysProxy.hpp"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <stack>
#include <tuple>

namespace cldes {

namespace op {
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
DESystem<NEvents, StorageIndex>
Synchronize(DESystemBase<NEvents, StorageIndex> const& aSys0,
            DESystemBase<NEvents, StorageIndex> const& aSys1)
{
    DESystem<NEvents, StorageIndex> sys = DESystem<NEvents, StorageIndex>(
      SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 });

    return sys;
}

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
SyncSysProxy<NEvents, StorageIndex>
SynchronizeStage1(DESystemBase<NEvents, StorageIndex> const& aSys0,
                  DESystemBase<NEvents, StorageIndex> const& aSys1)
{
    return SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 };
}

/*! \brief Transform a empty virtual system in a real system
 *
 * Calculate transitions and other parameters of a virtual system which
 * represents a virtual parallel composition.
 *
 * TODO: Cache lazy calculated transitions when the user visit transitions and
 * use it on synchronize stage 2
 *
 * @param aVirtualSys Reference to the system which will be transformed.
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 */
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeEmptyStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    // Estimated sparcity pattern: Calculate on the same loop than statesmap
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_;

    // Reserve space for transitions
    aVirtualSys.ResizeStatesEvents(aVirtualSys.states_number_);
    aVirtualSys.triplet_.reserve(sparcitypattern);
    aVirtualSys.bittriplet_.reserve(sparcitypattern +
                                    aVirtualSys.states_number_);

    // Calculate transitions
    for (StorageIndex qfrom = 0; qfrom < aVirtualSys.states_number_; ++qfrom) {
        aVirtualSys.bittriplet_.push_back(BitTriplet(qfrom, qfrom, true));

        aVirtualSys.SetStateEvents(qfrom, aVirtualSys.GetStateEvents(qfrom));
        aVirtualSys.SetInvStateEvents(qfrom,
                                      aVirtualSys.GetInvStateEvents(qfrom));

        auto event = 0u;
        auto event_it = aVirtualSys.GetStateEvents(qfrom);
        while (event_it != 0) {
            if (event_it.test(0)) {
                auto const qto = aVirtualSys.Trans(qfrom, event);

                aVirtualSys.triplet_.push_back(Triplet<NEvents>(
                  qfrom, qto, EventsSet<NEvents>{ 1ul << event }));

                if (qfrom != static_cast<StorageIndex>(qto)) {
                    aVirtualSys.bittriplet_.push_back(
                      BitTriplet(qto, qfrom, true));
                }
            }
            ++event;
            event_it >>= 1;
        }
    }

    return;
}

/*! \brief Transform a virtual system in a real system: optmized to supervisor
 * synthesis
 *
 * Calculate transitions and other parameters of a virtual system which
 * represents a virtual parallel composition. It does not calculate the
 * states_events_ table, since supervisors are final system and may have a big
 * size
 *
 * @param aVirtualSys Reference to the system which will be transformed.
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 */
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    SparseStatesMap<StorageIndex> statesmap;

    aVirtualSys.SetStatesNumber(aVirtualSys.virtual_states_.size());

    // Reserve space for transitions
    aVirtualSys.bittriplet_.reserve(aVirtualSys.states_number_);

    // Estimated sparcity pattern: Calculate on the same loop than statesmap
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_;

    // Map states to its new index
    // TODO: It SHOULD be returned in the future in a pair
    StorageIndex cst = 0;
    for (StorageIndex s : aVirtualSys.virtual_states_) {
        statesmap[s] = cst;
        aVirtualSys.bittriplet_.push_back(BitTriplet(cst, cst, true));
        ++cst;
    }

    // virtual_states_ is not necessary anymore and it can be a large vector
    aVirtualSys.virtual_states_.clear();

    // Remap marked states
    for (StorageIndex s0 : aVirtualSys.sys0_.GetMarkedStates()) {
        for (StorageIndex s1 : aVirtualSys.sys1_.GetMarkedStates()) {
            StorageIndex const key = s1 * aVirtualSys.n_states_sys0_ + s0;
            if (statesmap.contains(key)) {
                aVirtualSys.InsertMarkedState(statesmap[key]);
            }
        }
    }

    // Reserve space for transitions
    aVirtualSys.triplet_.reserve(sparcitypattern);
    aVirtualSys.bittriplet_.reserve(sparcitypattern +
                                    aVirtualSys.states_number_);

    // Calculate transitions
    while (!aVirtualSys.transtriplet_.empty()) {
        auto q_trans = aVirtualSys.transtriplet_.back();
        auto const q = q_trans.first;

        while (!q_trans.second->empty()) {
            auto const qto_e = q_trans.second->back();
            auto const qto = qto_e.first;

            if (statesmap.contains(qto)) {
                auto const event = qto_e.second;
                auto const qto_mapped = statesmap[qto];
                auto const q_mapped = statesmap[q];

                aVirtualSys.triplet_.push_back(Triplet<NEvents>(
                  q_mapped, qto_mapped, EventsSet<NEvents>{ 1ul << event }));
                aVirtualSys.bittriplet_.push_back(
                  BitTriplet(qto_mapped, q_mapped, true));
            }
            q_trans.second->pop_back();
        }
        delete q_trans.second;
        aVirtualSys.transtriplet_.pop_back();
    }

    return;
}

/*! \brief Remove bad states recursively
 *
 * Remove a state and all the states that become a bad state when the previous
 * one is removed.
 *
 * @param aVirtualSys reference to the virtual system which will have a state
 * removed
 * @param aC A hash table containing the states currently added to the
 * virtual system.
 * @param aQ The state to be removed.
 * @param aNonContrBit Bit array containing the virtual system non controllable
 * events.
 * @param aRmTable A hash table containing all the removed states so far
 */
template<uint8_t NEvents, typename StorageIndex>
void
RemoveBadStates(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys,
                TransMap<StorageIndex>& aC,
                StorageIndex const& aQ,
                EventsSet<NEvents> const& aNonContrBit,
                StatesTableHost<StorageIndex>& aRmTable)
{
    StatesStack<StorageIndex> f;
    f.push(aQ);
    aRmTable.insert(aQ);

    while (!f.empty()) {
        auto const x = f.top();
        f.pop();

        auto q_events = aVirtualSys.GetInvStateEvents(x);

        q_events &= aNonContrBit;

        cldes::ScalarType event = 0;
        while (q_events.any()) {
            if (q_events.test(0)) {
                StatesArray<StorageIndex> const finv =
                  aVirtualSys.InvTrans(x, event);

                for (StorageIndex s : finv) {
                    if (!aRmTable.contains(s)) {
                        f.push(s);
                        aRmTable.insert(s);
                        if (aC.contains(s)) {
                            delete aC[s];
                            aC.erase(s);
                        }
                    }
                }
            }
            ++event;
            q_events >>= 1;
        }
    }
    return;
}

/*! \brief Returns the monolithic supervisor of a plant and a pec
 *
 * @param aP Plant system const reference
 * @param aE Specs system const reference
 * @param aNonContr Hash table containing all non-controllable events indexes.
 */
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(DESystemBase<NEvents, StorageIndex> const& aP,
                DESystemBase<NEvents, StorageIndex> const& aE,
                EventsTableHost const& aNonContr)
{
    // Define new systems params: Stage1 is not necessary
    SyncSysProxy<NEvents, StorageIndex> virtualsys{ aP, aE };

    // non_contr in a bitarray structure
    EventsSet<NEvents> non_contr_bit;
    EventsSet<NEvents> p_non_contr_bit;

    // Evaluate which non contr event is in system and convert it to a
    // bitarray
    for (cldes::ScalarType event : aNonContr) {
        if (aP.GetEvents().test(event)) {
            p_non_contr_bit.set(event);
            if (virtualsys.events_.test(event)) {
                non_contr_bit.set(event);
            }
        }
    }

    // Supervisor states
    TransMap<StorageIndex> c;
    StatesTableHost<StorageIndex> rmtable;

    // f is a stack of states accessed in a dfs
    StatesStack<StorageIndex> f;

    // Initialize f and ftable with the initial state
    f.push(virtualsys.init_state_);

    // Allocate inverted graph, since we are search for inverse transitions
    virtualsys.AllocateInvertedGraph();

    while (!f.empty()) {
        auto const q = f.top();
        f.pop();

        if (!rmtable.contains(q) && !c.contains(q)) {
            auto const qx = q % virtualsys.n_states_sys0_;
            auto const q_events = virtualsys.GetStateEvents(q);

            auto const in_ncqx = p_non_contr_bit & aP.GetStateEvents(qx);
            auto const in_ncqx_and_q = in_ncqx & q_events;

            if (in_ncqx_and_q != in_ncqx) {
                // TODO: Fix template implicit instantiation
                RemoveBadStates<NEvents, StorageIndex>(
                  virtualsys, c, q, non_contr_bit, rmtable);
            } else {
                c[q] = new InvArgTrans<StorageIndex>();

                cldes::ScalarType event = 0;
                auto event_it = q_events;
                while (event_it.any()) {
                    if (event_it.test(0)) {
                        auto const fsqe = virtualsys.Trans(q, event);

                        if (!rmtable.contains(fsqe)) {
                            if (!c.contains(fsqe)) {
                                f.push(fsqe);
                            }
                        }
                        c[q]->push_back(std::make_pair(fsqe, event));
                    }
                    ++event;
                    event_it >>= 1;
                }
            }
        }
    }

    rmtable.clear();

    virtualsys.ClearInvertedGraph();

    // Swap new system states and sort it
    virtualsys.virtual_states_.reserve(c.size());
    virtualsys.transtriplet_.reserve(c.size());
    for (auto tr : c) {
        virtualsys.virtual_states_.push_back(tr.first);
        virtualsys.transtriplet_.push_back(tr);
    }
    c.clear();

    // Finish synching without removed states
    // SynchronizeStage2(virtualsys);

    // Transform virtual sys in a real sys by forcing conversion
    auto sys = DESystem<NEvents, StorageIndex>(virtualsys);

    // Remove non-accessible and non-coaccessible states
    sys.Trim();

    // bye
    return sys;
}

/*! \brief Returns the monolithic supervisor of a plant and a pec
 *
 * Calculate it in a lazy way.
// Assumes that vectors have length >= 1
 *
 * @param aP Plant system const reference
 * @param aE Specs system const reference
 * @param aNonContr Hash table containing all non-controllable events indexes.
 */
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(DESVector<NEvents, StorageIndex> const& aPlants,
                DESVector<NEvents, StorageIndex> const& aSpecs,
                EventsTableHost const& aNonContr)
{
    using SyncSysProxy = SyncSysProxy<NEvents, StorageIndex>;
    using DESystemBase = DESystemBase<NEvents, StorageIndex>;
    using DESystem = DESystem<NEvents, StorageIndex>;

    std::shared_ptr<DESystemBase const> plant[aPlants.size()];
    plant[0] = std::make_shared<DESystem const>(aPlants[0]);

    if (aPlants.size() > 1) {
        for (size_t k = 1; k < aPlants.size(); ++k) {
            std::shared_ptr<DESystemBase const> sync =
              std::make_shared<SyncSysProxy>(
                SyncSysProxy{ *(plant[k - 1]), aPlants[k] });
            plant[k] = sync;
        }
    }

    std::shared_ptr<DESystemBase const> spec[aPlants.size()];
    spec[0] = std::make_shared<DESystem>(aSpecs[0]);

    if (aSpecs.size() > 1) {
        for (size_t k = 1; k < aSpecs.size(); ++k) {
            std::shared_ptr<DESystemBase const> sync =
              std::make_shared<SyncSysProxy>(
                SyncSysProxy{ *(spec[k - 1]), aSpecs[k] });
            spec[k] = sync;
        }
    }

    auto const supervisor = SupervisorSynth<NEvents, StorageIndex>(
      *(plant[aPlants.size() - 1]), *(spec[aSpecs.size() - 1]), aNonContr);

    return supervisor;
}

} // namespace op
} // namespace cldes

#endif // DESYSTEM_HPP
