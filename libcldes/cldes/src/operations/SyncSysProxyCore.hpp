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
 Description: Definition of SyncSysProxy class.
 =========================================================================
*/

#include <algorithm>

template<uint8_t NEvents, typename StorageIndex>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::SyncSysProxy(
  cldes::DESystemBase<NEvents, StorageIndex> const& aSys0,
  cldes::DESystemBase<NEvents, StorageIndex> const& aSys1)
  : sys0_{ aSys0 }
  , sys1_{ aSys1 }
{
    n_states_sys0_ = aSys0.GetStatesNumber();
    this->states_number_ = n_states_sys0_ * aSys1.GetStatesNumber();
    this->init_state_ = aSys1.GetInitialState() * aSys0.GetStatesNumber() +
                        aSys0.GetInitialState();

    auto const in_both = aSys0.GetEvents() & aSys1.GetEvents();

    only_in_0_ = aSys0.GetEvents() ^ in_both;
    only_in_1_ = aSys1.GetEvents() ^ in_both;

    // Initialize events_ from base class
    this->events_ = aSys0.GetEvents() | aSys1.GetEvents();

    for (auto q0 : aSys0.GetMarkedStates()) {
        for (auto q1 : aSys1.GetMarkedStates()) {
            this->marked_states_.emplace(q1 * n_states_sys0_ + q0);
        }
    }
}

template<uint8_t NEvents, typename StorageIndex>
bool
cldes::op::SyncSysProxy<NEvents, StorageIndex>::ContainsTrans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsTrans(qy, aEvent);

    auto contains = false;

    if ((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
        (in_y && only_in_1_.test(aEvent))) {
        contains = true;
    }

    return contains;
}

template<uint8_t NEvents, typename StorageIndex>
typename cldes::op::SyncSysProxy<NEvents, StorageIndex>::StorageIndexSigned
cldes::op::SyncSysProxy<NEvents, StorageIndex>::Trans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsTrans(qy, aEvent);

    if (!((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
          (in_y && only_in_1_.test(aEvent)))) {
        return -1;
    }

    if (in_x && in_y) {
        auto const q0 = sys0_.Trans(qx, aEvent);
        auto const q1 = sys1_.Trans(qy, aEvent);

        return q1 * n_states_sys0_ + q0;
    } else if (in_x) {
        auto const q0 = sys0_.Trans(qx, aEvent);
        return qy * n_states_sys0_ + q0;
    } else { // in_y
        auto const q1 = sys1_.Trans(qy, aEvent);
        return q1 * n_states_sys0_ + qx;
    }
}

template<uint8_t NEvents, typename StorageIndex>
bool
cldes::op::SyncSysProxy<NEvents, StorageIndex>::ContainsInvTrans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

    auto contains = false;

    if ((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
        (in_y && only_in_1_.test(aEvent))) {
        contains = true;
    }

    return contains;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::StatesArray<StorageIndex>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::InvTrans(
  StorageIndex const& aQ,
  ScalarType const& aEvent) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

    StatesArray<StorageIndex> inv_transitions;

    if (!((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
          (in_y && only_in_1_.test(aEvent)))) {
        return inv_transitions;
    }

    if (in_x && in_y) {
        auto const inv_trans_0 = sys0_.InvTrans(qx, aEvent);
        auto const inv_trans_1 = sys1_.InvTrans(qy, aEvent);

        inv_transitions.reserve(inv_trans_0.size() + inv_trans_1.size());

        for (auto q0 : inv_trans_0) {
            for (auto q1 : inv_trans_1) {
                auto const q_from = q1 * n_states_sys0_ + q0;
                inv_transitions.push_back(q_from);
            }
        }
    } else if (in_x) {
        auto const inv_trans_0 = sys0_.InvTrans(qx, aEvent);

        inv_transitions.reserve(inv_trans_0.size());

        for (auto q : inv_trans_0) {
            auto const q_from = qy * n_states_sys0_ + q;
            inv_transitions.push_back(q_from);
        }
    } else { // in_y
        auto const inv_trans_1 = sys1_.InvTrans(qy, aEvent);

        inv_transitions.reserve(inv_trans_1.size());

        for (auto q : inv_trans_1) {
            auto const q_from = q * n_states_sys0_ + qx;
            inv_transitions.push_back(q_from);
        }
    }

    return inv_transitions;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::
operator DESystem<NEvents, StorageIndex>()
{
    std::sort(virtual_states_.begin(), virtual_states_.end());
    SynchronizeStage2(*this);

    DESystem<NEvents, StorageIndex> sys{};
    sys.states_number_ = this->states_number_;

    // Resize adj matrices
    sys.graph_.resize(this->states_number_, this->states_number_);
    sys.bit_graph_.resize(this->states_number_, this->states_number_);

    // Move triplets to graph storage
    sys.graph_.setFromTriplets(triplet_.begin(), triplet_.end());
    sys.bit_graph_.setFromTriplets(
      bittriplet_.begin(), bittriplet_.end(), [](bool const&, bool const&) {
          return true;
      });

    triplet_.clear();
    bittriplet_.clear();

    sys.graph_.makeCompressed();
    sys.bit_graph_.makeCompressed();

    return sys;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::EventsSet<NEvents>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::GetStateEvents(
  StorageIndex const& aQ) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const state_event_0 = sys0_.GetStateEvents(qx);
    auto const state_event_1 = sys1_.GetStateEvents(qy);
    auto const state_event = (state_event_0 & state_event_1) |
                             (state_event_0 & only_in_0_) |
                             (state_event_1 & only_in_1_);

    return state_event;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::EventsSet<NEvents>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::GetInvStateEvents(
  StorageIndex const& aQ) const
{
    // q = (qx, qy)
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const state_event_0 = sys0_.GetInvStateEvents(qx);
    auto const state_event_1 = sys1_.GetInvStateEvents(qy);
    auto const state_event = (state_event_0 & state_event_1) |
                             (state_event_0 & only_in_0_) |
                             (state_event_1 & only_in_1_);

    return state_event;
}
