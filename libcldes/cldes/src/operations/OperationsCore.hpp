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
 Description: Definition of operation functions: Parallel composition,
 virtual parallel composition and supervisor synthesis.
 =========================================================================
*/

#include "cldes/DESystem.hpp"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>

template<uint8_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>
cldes::op::Synchronize(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                       cldes::DESystem<NEvents, StorageIndex> const& aSys1)
{
    DESystem<NEvents, StorageIndex> sys = DESystem<NEvents, StorageIndex>(
      SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 });

    return sys;
}

template<uint8_t NEvents, typename StorageIndex>
typename cldes::op::SyncSysProxy<NEvents, StorageIndex>
cldes::op::SynchronizeStage1(
  cldes::DESystem<NEvents, StorageIndex> const& aSys0,
  cldes::DESystem<NEvents, StorageIndex> const& aSys1)
{
    return SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 };
}

template<uint8_t NEvents, typename StorageIndex>
void
cldes::op::SynchronizeEmptyStage2(
  cldes::op::SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
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

template<uint8_t NEvents, typename StorageIndex>
void
cldes::op::SynchronizeStage2(
  cldes::op::SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    SparseStatesMap<StorageIndex> statesmap;

    aVirtualSys.SetStatesNumber(aVirtualSys.virtual_states_.size());

    // Reserve space for transitions
    aVirtualSys.bittriplet_.reserve(aVirtualSys.states_number_);

    // Estimated sparcity pattern: Calculate on the same loop than statesmap
    StorageIndex sparcitypattern = 0;

    // Map states to its new index
    // TODO: It SHOULD be returned in the future in a pair
    StorageIndex cst = 0;
    for (StorageIndex s : aVirtualSys.virtual_states_) {
        statesmap[s] = cst;
        sparcitypattern += aVirtualSys.GetStateEvents(s).count();
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

template<uint8_t NEvents, typename StorageIndex>
void
cldes::op::RemoveBadStates(
  cldes::op::SyncSysProxy<NEvents, StorageIndex>& aVirtualSys,
  TransMap<StorageIndex>& aC,
  StorageIndex const& aQ,
  EventsSet<NEvents> const& aNonContrBit,
  cldes::op::StatesTableHost<StorageIndex>& aRmTable)
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

template<uint8_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>
cldes::op::SupervisorSynth(cldes::DESystem<NEvents, StorageIndex> const& aP,
                           cldes::DESystem<NEvents, StorageIndex> const& aE,
                           op::EventsTableHost const& aNonContr)
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
