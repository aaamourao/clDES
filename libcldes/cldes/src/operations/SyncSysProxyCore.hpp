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

template<uint8_t NEvents, typename StorageIndex>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::SyncSysProxy(
  cldes::DESystemBase<NEvents, StorageIndex> const& aSys0,
  cldes::DESystemBase<NEvents, StorageIndex> const& aSys1)
  : DESystemBase<NEvents, StorageIndex>{
      aSys0.states_number_ * aSys1.states_number_,
      aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_
  }
{
    sys0_ = aSys0;
    sys1_ = aSys1;
    n_states_sys0_ = aSys0.states_number_;

    auto const in_both = aSys0.events_ & aSys1.events_;

    only_in_0_ = aSys0.events_ ^ in_both;
    only_in_1_ = aSys1.events_ ^ in_both;

    // Initialize events_ from base class
    this->events_ = aSys0.events_ & aSys1.events_;
}

template<uint8_t NEvents, typename StorageIndex>
bool
cldes::op::SyncSysProxy<NEvents, StorageIndex>::ContainsTrans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent) const
{
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsTrans(qy, aEvent);

    auto contains = false;

    if ((in_x && in_y) || (in_x && only_in_0_) || (in_y && only_in_1_)) {
        contains = true;
    }

    return contains;
}

template<uint8_t NEvents, typename StorageIndex>
typename cldes::op::SyncSysProxy<NEvents, StorageIndex>::StorageIndexSigned
cldes::op::SyncSysProxy<NEvents, StorageIndex>::Trans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent)
{
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsTrans(qy, aEvent);

    if (!(in_x && in_y) || !(in_x && only_in_0_) || !(in_y && only_in_1_)) {
        return -1;
    }

    if (in_x && in_y) {
        auto const trans_0 = sys0_.Trans(qx, aEvent);
        auto const trans_1 = sys1_.Trans(qy, aEvent);

        for (auto q0 : trans_0) {
            for (auto q1 : trans_1) {
                return q1 * n_states_sys0_ + q0;
            }
        }
    } else if (in_x) {
        auto const trans_0 = sys0_.Trans(qx);

        for (auto q : trans_0) {
            return qy * n_states_sys0_ + q;
        }
    } else { // in_y
        auto const trans_1 = sys1_.Trans(qy);

        for (auto q : trans_1) {
            return q * n_states_sys0_ + qx;
        }
    }
}

template<uint8_t NEvents, typename StorageIndex>
bool
cldes::op::SyncSysProxy<NEvents, StorageIndex>::ContainsInvTrans(
  StorageIndex const& aQ,
  cldes::ScalarType const& aEvent) const
{
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

    auto contains = false;

    if ((in_x && in_y) || (in_x && only_in_0_) || (in_y && only_in_1_)) {
        contains = true;
    }

    return contains;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::StatesArray<StorageIndex>
cldes::op::SyncSysProxy<NEvents, StorageIndex>::InvTrans(
  StorageIndex const& aQ,
  ScalarType const& aEvent)
{
    auto const qx = aQ % n_states_sys0_;
    auto const qy = aQ / n_states_sys0_;

    auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
    auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

    StatesArray<StorageIndex> inv_transitions;

    if (!(in_x && in_y) || !(in_x && only_in_0_) || !(in_y && only_in_1_)) {
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
        auto const inv_trans_0 = sys0_.InvTrans(qx);

        inv_transitions.reserve(inv_trans_0.size());

        for (auto q : inv_trans_0) {
            auto const q_from = qy * n_states_sys0_ + q;
            inv_transitions.push_back(q_from);
        }
    } else { // in_y
        auto const inv_trans_1 = sys1_.InvTrans(qy);

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
    return DESystem<NEvents, StorageIndex>{};
}
