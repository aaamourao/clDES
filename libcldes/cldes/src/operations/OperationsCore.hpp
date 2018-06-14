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
    using Triplet = cldes::Triplet<NEvents>;
    using RowIterator = Eigen::InnerIterator<
      typename DESystem<NEvents, StorageIndex>::GraphHostData const>;

    auto const in_both = aSys0.events_ & aSys1.events_;
    auto const only_in_0 = aSys0.events_ ^ in_both;
    auto const only_in_1 = aSys1.events_ ^ in_both;

    // Calculate new marked states
    typename DESystem<NEvents, StorageIndex>::StatesSet marked_states;
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            marked_states.insert(s1 * aSys0.states_number_ + s0);
        }
    }

    // Create new system without transitions
    DESystem<NEvents, StorageIndex> sys{
        aSys0.states_number_ * aSys1.states_number_,
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_,
        marked_states
    };

    // Set private params
    sys.events_ = aSys0.events_ | aSys1.events_;

    // Alias to states_number_
    StorageIndex const nstates = sys.states_number_;

    // Calculate sparcity pattern
    StorageIndex const sparcitypattern = sys.events_.count() * nstates;

    // Reserve space for transitions
    std::vector<Triplet> triplet;
    std::vector<BitTriplet> bittriplet;

    triplet.reserve(sparcitypattern);
    bittriplet.reserve(sparcitypattern);

    // Calculate transitions
    for (StorageIndex q = 0; q < nstates; ++q) {
        auto const qx = q % aSys0.states_number_;
        auto const qy = q / aSys0.states_number_;

        // Calculate sys inverse states events
        sys.inv_states_events_[q] =
          (aSys0.inv_states_events_[qx] & aSys1.inv_states_events_[qy]) |
          (aSys0.inv_states_events_[qx] & only_in_0) |
          (aSys1.inv_states_events_[qy] & only_in_1);

        // Calculate sys states events
        auto q_events = (aSys0.states_events_[qx] & aSys1.states_events_[qy]) |
                        (aSys0.states_events_[qx] & only_in_0) |
                        (aSys1.states_events_[qy] & only_in_1);
        sys.states_events_[q] = q_events;

        // Add loop to bit_graph_ : bit graph = graph.in_bits + identity
        bittriplet.push_back(BitTriplet(q, q, true));

        cldes::ScalarType event = 0;
        while (q_events.any()) {
            if (q_events.test(0)) {
                StorageIndex qto;

                StorageIndex xto = 0;
                StorageIndex yto = 0;

                auto const is_in_p = aSys0.events_.test(event);
                auto const is_in_e = aSys1.events_.test(event);

                if (is_in_p && is_in_e) {
                    for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
                        if (pe.value().test(event)) {
                            xto = pe.col();
                            break;
                        }
                    }
                    for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
                        if (ee.value().test(event)) {
                            yto = ee.col();
                            break;
                        }
                    }
                } else if (is_in_e) {
                    for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
                        if (ee.value().test(event)) {
                            xto = qx;
                            yto = ee.col();
                            break;
                        }
                    }
                } else {
                    for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
                        if (pe.value().test(event)) {
                            xto = pe.col();
                            yto = qy;
                            break;
                        }
                    }
                }

                qto = yto * aSys0.states_number_ + xto;

                triplet.push_back(
                  Triplet(q, qto, EventsSet<NEvents>{ 1ul << event }));
                if (q != qto) {
                    bittriplet.push_back(BitTriplet(qto, q, true));
                }
            }
            ++event;
            q_events >>= 1;
        }
    }

    // Remove aditional space
    sys.graph_.setFromTriplets(triplet.begin(), triplet.end());
    sys.bit_graph_.setFromTriplets(
      bittriplet.begin(), bittriplet.end(), [](bool const&, bool const&) {
          return true;
      });

    sys.graph_.makeCompressed();
    sys.bit_graph_.makeCompressed();

    return sys;
}

template<uint8_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>
cldes::op::SynchronizeStage1(
  cldes::DESystem<NEvents, StorageIndex> const& aSys0,
  cldes::DESystem<NEvents, StorageIndex> const& aSys1)
{
    auto const in_both = aSys0.events_ & aSys1.events_;
    auto const only_in_0 = aSys0.events_ ^ in_both;
    auto const only_in_1 = aSys1.events_ ^ in_both;

    // Calculate new marked states
    typename DESystem<NEvents, StorageIndex>::StatesSet marked_states;
    for (StorageIndex s0 : aSys0.marked_states_) {
        for (StorageIndex s1 : aSys1.marked_states_) {
            marked_states.insert(s1 * aSys0.states_number_ + s0);
        }
    }

    // Create new system without transitions
    DESystem<NEvents, StorageIndex> virtualsys{
        aSys0.states_number_ * aSys1.states_number_,
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_,
        marked_states
    };

    // New system params
    virtualsys.states_events_.reserve(aSys0.states_number_ *
                                      aSys1.states_number_);
    virtualsys.inv_states_events_.reserve(aSys0.states_number_ *
                                          aSys1.states_number_);

    // Calculate params
    for (StorageIndex ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
        for (StorageIndex ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
            auto const key = ix1 * aSys0.states_number_ + ix0;

            virtualsys.virtual_states_.push_back(key);

            virtualsys.states_events_[key] =
              (aSys0.states_events_[ix0] & aSys1.states_events_[ix1]) |
              (aSys0.states_events_[ix0] & only_in_0) |
              (aSys1.states_events_[ix1] & only_in_1);
            virtualsys.inv_states_events_[key] =
              (aSys0.inv_states_events_[ix0] & aSys1.inv_states_events_[ix1]) |
              (aSys0.inv_states_events_[ix0] & only_in_0) |
              (aSys1.inv_states_events_[ix1] & only_in_1);
        }
    }

    // Set private params
    virtualsys.events_ = aSys0.events_ | aSys1.events_;

    return virtualsys;
}

template<uint8_t NEvents, typename StorageIndex>
void
cldes::op::SynchronizeStage2(
  cldes::op::SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    SparseStatesMap<StorageIndex> statesmap;

    aVirtualSys.SetStatesNumber(aVirtualSys.virtual_states_.size());

    // Estimate sparcity pattern
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_;

    // Reserve space for transitions
    aVirtualSys.bittriplet_.reserve(sparcitypattern +
                                    aVirtualSys.states_number_);

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
