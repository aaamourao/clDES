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

 File: operations.cpp
 Description: Definition of operation functions.
 =========================================================================
*/

#include <algorithm>
#include <cmath>
// #include "backend/oclbackend.hpp"
#include "desystem.hpp"
// #include "desystemcl.hpp"
#include "transition_proxy.hpp"
#include <Eigen/Sparse>
#include <QtCore/QVector>

/*
DESystemCL op::Synchronize(DESystemCL &aSys0, DESystemCL &aSys1) {
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
        CL_MEM_READ_WRITE, table_size * sizeof(StatesTuple), nullptr);

    auto syncstage1kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage1");

    // Set Work groups size
    SetWorkGroups_(&syncstage1kernel, aSys0.states_number_,
                   aSys1.states_number_, 1, 1);

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage1kernel(
        states_tuple_dev, static_cast<cl_uint>(aSys0.states_number_)));

    // Get the result and saves it on host memory
    auto states_tuple_host = new StatesTuple[table_size];
    clEnqueueReadBuffer(aSys0.backend_ptr_->CommandQueue(), states_tuple_dev,
                        CL_TRUE, 0, table_size * sizeof(StatesTuple),
                        (void *)states_tuple_host, 0, NULL, NULL);

    if (aSys0.is_cache_outdated_) {
        aSys0.UpdateGraphCache_();
    }

    if (aSys1.is_cache_outdated_) {
        aSys1.UpdateGraphCache_();
    }

    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    auto syncstage2kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage2");

    // Set Work groups size
    SetWorkGroups_(&syncstage2kernel, table_size, aSys0.events_.size(), 1, 1);

    auto asys0_events = CalcEventsInt_(aSys0.events_);
    auto asys1_events = CalcEventsInt_(aSys1.events_);

    auto gcd_private = CalcGCD_(asys0_events, asys1_events);
    auto asys0_private = asys0_events / gcd_private;
    auto asys1_private = asys1_events / gcd_private;

    viennacl::matrix<float> result_dev(table_size, table_size);
    result_dev.clear();

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage2kernel(
        states_tuple_dev, aSys0.device_graph_->handle1().opencl_handle(),
        aSys0.device_graph_->handle2().opencl_handle(),
        aSys0.device_graph_->handle().opencl_handle(),
        static_cast<cl_uint>(aSys0.device_graph_->rows()), asys0_private,
        aSys1.device_graph_->handle1().opencl_handle(),
        aSys1.device_graph_->handle2().opencl_handle(),
        aSys1.device_graph_->handle().opencl_handle(), asys1_private,
        result_dev.handle().opencl_handle(),
        static_cast<cl_uint>(result_dev.internal_rows())));

    // Copy device graph to host memory
    DESystemCL sync_sys(table_size, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));

    return sync_sys;
}
*/

cldes::DESystem
cldes::op::Synchronize(cldes::DESystem const& aSys0,
                       cldes::DESystem const& aSys1)
{
    using Triplet = Eigen::Triplet<EventsBitArray>;
    using BitTriplet = Eigen::Triplet<bool>;
    using RowIterator = Eigen::InnerIterator<DESystem::GraphHostData const>;

    auto const in_both = aSys0.events_ & aSys1.events_;
    auto const only_in_0 = aSys0.events_ ^ in_both;
    auto const only_in_1 = aSys1.events_ ^ in_both;

    // Calculate new marked states
    DESystem::StatesSet marked_states;
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            marked_states.insert(s1 * aSys0.states_number_ + s0);
        }
    }

    // Create new system without transitions
    DESystem sys{ aSys0.states_number_ * aSys1.states_number_,
                  aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_,
                  marked_states };

    // Set private params
    sys.events_ = aSys0.events_ | aSys1.events_;

    // Alias to states_number_
    auto const nstates = sys.states_number_;

    // Calculate sparcity pattern
    auto const sparcitypattern = sys.events_.count() * nstates;

    // Reserve space for transitions
    std::vector<Triplet> triplet;
    std::vector<BitTriplet> bittriplet;

    triplet.reserve(sparcitypattern);
    bittriplet.reserve(sparcitypattern);

    // Calculate transitions
    for (auto q = 0ul; q < nstates; ++q) {
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

        auto event = 0ul;
        while (q_events.any()) {
            if (q_events.test(0)) {
                int qto;

                int xto = -1;
                int yto = -1;

                bool const is_in_p = aSys0.events_.test(event);
                bool const is_in_e = aSys1.events_.test(event);

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
                  Triplet(q, qto, EventsBitArray{ 1ul << event }));
                if (q != static_cast<unsigned long>(qto)) {
                    bittriplet.push_back(BitTriplet(qto, q, true));
                }
            }
            ++event;
            q_events >>= 1ul;
        }
    }

    // Remove aditional space
    sys.graph_.setFromTriplets(triplet.begin(), triplet.end());
    sys.bit_graph_.setFromTriplets(bittriplet.begin(), bittriplet.end());

    sys.graph_.makeCompressed();
    sys.bit_graph_.makeCompressed();

    return sys;
}

/*
op::StatesTable *op::SynchronizeStage1(DESystemCL const &aSys0,
                                       DESystemCL const &aSys1) {
    SynchronizeStage2(syncsys, aSys0, aSys1);
    return syncsys;
}
*/

cldes::DESystem
cldes::op::SynchronizeStage1(cldes::DESystem const& aSys0,
                             cldes::DESystem const& aSys1)
{
    auto const in_both = aSys0.events_ & aSys1.events_;
    auto const only_in_0 = aSys0.events_ ^ in_both;
    auto const only_in_1 = aSys1.events_ ^ in_both;

    // Calculate new marked states
    DESystem::StatesSet marked_states;
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            marked_states.insert(s1 * aSys0.states_number_ + s0);
        }
    }

    // Create new system without transitions
    DESystem virtualsys{ aSys0.states_number_ * aSys1.states_number_,
                         aSys1.init_state_ * aSys0.states_number_ +
                           aSys0.init_state_,
                         marked_states };

    // New system params
    virtualsys.states_events_.reserve(aSys0.states_number_ *
                                      aSys1.states_number_);
    virtualsys.inv_states_events_.reserve(aSys0.states_number_ *
                                          aSys1.states_number_);

    // Calculate params
    for (auto ix0 = 0ul; ix0 < aSys0.states_number_; ++ix0) {
        for (auto ix1 = 0ul; ix1 < aSys1.states_number_; ++ix1) {
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

/*
op::StatesTable *op::SynchronizeStage1(DESystemCL const &aSys0,
                                       DESystemCL const &aSys1) {
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
        CL_MEM_WRITE_ONLY, table_size * sizeof(StatesTuple), nullptr);

    auto syncstage1kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage1");

    // Set Work groups size
    SetWorkGroups_(&syncstage1kernel, aSys0.states_number_,
                   aSys1.states_number_, 1, 1);

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage1kernel(
        states_tuple_dev, static_cast<cl_uint>(aSys0.states_number_)));

    // Get the result and saves it on host memory
    auto states_tuple_host = new StatesTuple[table_size];
    clEnqueueReadBuffer(aSys0.backend_ptr_->CommandQueue(), states_tuple_dev,
                        CL_TRUE, 0, table_size * sizeof(StatesTuple),
                        (void *)states_tuple_host, 0, NULL, NULL);

    auto states_table = new StatesTable;
    states_table->tsize = table_size;
    states_table->table = states_tuple_host;

    return states_table;
}

DESystemCL op::SynchronizeStage2(op::StatesTable const *aTable,
                                 DESystemCL &aSys0, DESystemCL &aSys1) {
    if (aSys0.is_cache_outdated_) {
        aSys0.UpdateGraphCache_();
    }

    if (aSys1.is_cache_outdated_) {
        aSys1.UpdateGraphCache_();
    }

    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
        CL_MEM_READ_ONLY, aTable->tsize * sizeof(StatesTuple), aTable->table);

    auto syncstage2kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage2");

    // Set Work groups size
    SetWorkGroups_(&syncstage2kernel, aTable->tsize, aSys0.events_.size(), 1,
                   1);

    auto asys0_events = CalcEventsInt_(aSys0.events_);
    auto asys1_events = CalcEventsInt_(aSys1.events_);

    auto gcd_private = CalcGCD_(asys0_events, asys1_events);
    auto asys0_private = asys0_events / gcd_private;
    auto asys1_private = asys1_events / gcd_private;

    viennacl::matrix<float> result_dev(aTable->tsize, aTable->tsize);
    result_dev.clear();

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage2kernel(
        states_tuple_dev, aSys0.device_graph_->handle1().opencl_handle(),
        aSys0.device_graph_->handle2().opencl_handle(),
        aSys0.device_graph_->handle().opencl_handle(),
        static_cast<cl_uint>(aSys0.device_graph_->rows()), asys0_private,
        aSys1.device_graph_->handle1().opencl_handle(),
        aSys1.device_graph_->handle2().opencl_handle(),
        aSys1.device_graph_->handle().opencl_handle(), asys1_private,
        result_dev.handle().opencl_handle(),
        static_cast<cl_uint>(result_dev.internal_rows())));

    // Copy device graph to host memory
    DESystemCL sync_sys(aTable->tsize, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));

    return sync_sys;
}
*/

void
cldes::op::SynchronizeStage2(cldes::DESystem& aVirtualSys,
                             cldes::DESystem const& aSys0,
                             cldes::DESystem const& aSys1)
{
    using Triplet = Eigen::Triplet<EventsBitArray>;
    using BitTriplet = Eigen::Triplet<bool>;

    // Alias to new size
    auto const nstates = aVirtualSys.virtual_states_.size();

    // Update new size
    aVirtualSys.states_number_ = nstates;

    aVirtualSys.events_ = aSys0.events_ | aSys1.events_;

    // Resize adj matrices if necessary
    aVirtualSys.states_events_.reserve(nstates);
    aVirtualSys.inv_states_events_.reserve(nstates);
    aVirtualSys.graph_.resize(nstates, nstates);
    aVirtualSys.bit_graph_.resize(nstates, nstates);

    // Estimate sparcity pattern
    auto const sparcitypattern = aVirtualSys.events_.count() * nstates;

    // Reserve space for transitions
    std::vector<Triplet> triplet;
    std::vector<BitTriplet> bittriplet;

    triplet.reserve(sparcitypattern);
    bittriplet.reserve(sparcitypattern + aVirtualSys.states_number_);

    QHash<cldes_size_t, cldes_size_t> statesmap;
    statesmap.reserve(aVirtualSys.virtual_states_.size());
    auto cst = 0ul;
    foreach (cldes_size_t s, aVirtualSys.virtual_states_) {
        statesmap[s] = cst;
        bittriplet.push_back(BitTriplet(cst, cst, true));

        ++cst;
    }

    // Calculate transitions
    while (!aVirtualSys.transtriplet_.empty()) {
        auto q_trans = aVirtualSys.transtriplet_.back();
        auto const q = q_trans.first;

        while (!q_trans.second.empty()) {
            auto const qto_e = q_trans.second.back();
            auto const qto = qto_e.first;

            if (statesmap.contains(qto)) {
                auto const event = qto_e.second;
                auto const qto_mapped = statesmap.value(qto);
                auto const q_mapped = statesmap.value(q);

                triplet.push_back(
                  Triplet(q_mapped, qto_mapped, EventsBitArray{ event }));
                bittriplet.push_back(BitTriplet(qto_mapped, q_mapped, true));
            }
            q_trans.second.pop_back();
        }
        aVirtualSys.transtriplet_.pop_back();
    }

    // Remove aditional space
    aVirtualSys.graph_.setFromTriplets(triplet.begin(), triplet.end());
    aVirtualSys.bit_graph_.setFromTriplets(bittriplet.begin(),
                                           bittriplet.end());
    aVirtualSys.graph_.makeCompressed();
    aVirtualSys.bit_graph_.makeCompressed();

    // Remap marked states
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            auto const key = s1 * aSys0.states_number_ + s0;
            if (statesmap.contains(key)) {
                aVirtualSys.marked_states_.insert(statesmap.value(key));
            }
        }
    }

    // It only works for init_state = 0;
    aVirtualSys.init_state_ =
      aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_;

    // Remove set that will not be used anymore
    aVirtualSys.virtual_states_.clear();
    aVirtualSys.only_in_0_.reset();
    aVirtualSys.only_in_1_.reset();
}

/*
bool op::ExistTransitionVirtual(DESystem const &aSys0, DESystem const
&aSys1, op::StatesTupleSTL const q, ScalarType const event) { bool const
is_in_p = aSys0.events_[event]; bool const is_in_e = aSys1.events_[event];

    bool const is_in_x = (aSys0.states_events_[q.first])[event];
    bool const is_in_y = (aSys1.states_events_[q.second])[event];

    bool exist_transition = false;

    if ((is_in_x && is_in_y) || (is_in_x && !is_in_e) ||
        (is_in_y && !is_in_p)) {
        exist_transition = true;
    }

    return exist_transition;
}
    */

cldes::cldes_size_t
cldes::op::TransitionVirtual(cldes::DESystem const& aSys0,
                             cldes::DESystem const& aSys1,
                             cldes::cldes_size_t const& q,
                             cldes::ScalarType const& event)
{
    using RowIterator = Eigen::InnerIterator<DESystem::GraphHostData const>;

    bool const is_in_p = aSys0.events_.test(event);
    bool const is_in_e = aSys1.events_.test(event);

    auto const qx = q % aSys0.states_number_;
    auto const qy = q / aSys0.states_number_;

    int xid = -1;
    int yid = -1;

    if (is_in_p && is_in_e) {
        for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                xid = pe.col();
                break;
            }
        }
        for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                yid = ee.col();
                break;
            }
        }
    } else if (is_in_e) {
        for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                xid = qx;
                yid = ee.col();
                break;
            }
        }
    } else {
        for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                xid = pe.col();
                yid = qy;
                break;
            }
        }
    }

    return yid * aSys0.states_number_ + xid;
}

// This function assumes that there is an inverse transition.
template<class EventsType>
static cldes::op::StatesArray
__TransitionVirtualInv(EventsType const& aEventsP,
                       EventsType const& aEventsE,
                       cldes::op::GraphType const& aInvGraphP,
                       cldes::op::GraphType const& aInvGraphE,
                       cldes::cldes_size_t const& q,
                       cldes::ScalarType const& event)
{
    using RowIterator =
      Eigen::InnerIterator<cldes::DESystem::GraphHostData const>;

    auto const qx = q % aInvGraphP.rows();
    auto const qy = q / aInvGraphP.rows();

    bool const is_in_p = aEventsP.test(event);
    bool const is_in_e = aEventsE.test(event);

    cldes::op::StatesArray ret;
    // ret.reserve(cldes::g_max_events);

    auto const p_size = aInvGraphP.rows();

    if (is_in_p && is_in_e) {
        cldes::op::StatesArray pstates;
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                pstates.push_back(pe.col());
            }
        }
        for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                foreach (cldes::cldes_size_t sp, pstates) {
                    ret.push_back(ee.col() * p_size + sp);
                }
            }
        }
    } else if (is_in_p) { // Is only in p: is_in_p && !is_in_e
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                ret.push_back(qy * p_size + pe.col());
            }
        }
    } else { // Is only in e: !is_in_p && is_in_e
        for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                ret.push_back(ee.col() * p_size + qx);
            }
        }
    }

    return ret;
}

void
cldes::op::RemoveBadStates(cldes::DESystem& aVirtualSys,
                           cldes::DESystem const& aP,
                           cldes::DESystem const& aE,
                           cldes::op::GraphType const& aInvGraphP,
                           cldes::op::GraphType const& aInvGraphE,
                           cldes::op::StatesTableSTL& C,
                           cldes::cldes_size_t const& q,
                           cldes::EventsBitArray const& bit_non_contr,
                           cldes::op::StatesTableSTL& rmtable)
{
    StatesStack f;
    f.push(q);
    rmtable.insert(q);

    while (!f.isEmpty()) {
        cldes_size_t const x = f.pop();

        auto const x0 = x % aInvGraphP.rows();
        auto const x1 = x / aInvGraphP.rows();

        auto q_events =
          (aP.inv_states_events_[x0] & aE.inv_states_events_[x1]) |
          (aP.inv_states_events_[x0] & aVirtualSys.only_in_0_) |
          (aE.inv_states_events_[x1] & aVirtualSys.only_in_1_);

        q_events &= bit_non_contr;

        auto event = 0ul;
        while (q_events.any()) {
            if (q_events.test(0)) {
                auto const finv = __TransitionVirtualInv(
                  aP.events_, aE.events_, aInvGraphP, aInvGraphE, x, event);

                foreach (cldes_size_t s, finv) {
                    if (!rmtable.contains(s)) {
                        f.push(s);
                        C.remove(s);
                        rmtable.insert(s);
                    }
                }
            }
            ++event;
            q_events >>= 1ul;
        }
    }
    return;
}

cldes::DESystem
cldes::op::SupervisorSynth(cldes::DESystem const& aP,
                           cldes::DESystem const& aE,
                           QSet<cldes::ScalarType> const& non_contr)
{
    DESystem::GraphHostData const p_invgraph = aP.graph_.transpose();
    DESystem::GraphHostData const e_invgraph = aE.graph_.transpose();

    // TODO: It may be in the new SynchronizeStage1
    DESystem virtualsys;
    virtualsys.init_state_ =
      aE.init_state_ * aP.states_number_ + aP.init_state_;
    virtualsys.is_cache_outdated_ = true;
    virtualsys.events_ = aP.events_ | aE.events_;
    // TODO

    // Alias to events in both systems
    auto const in_both = aP.events_ & aE.events_;

    // Calculate event parameters
    virtualsys.only_in_0_ = aP.events_ ^ in_both;
    virtualsys.only_in_1_ = aE.events_ ^ in_both;

    // non_contr in a bitarray structure
    EventsBitArray non_contr_bit;
    EventsBitArray p_non_contr_bit;

    // Evaluate which non contr event is in system and convert it to a
    // bitarray
    foreach (auto event, non_contr) {
        if (aP.events_.test(event)) {
            p_non_contr_bit.set(event);
            if (virtualsys.events_.test(event)) {
                non_contr_bit.set(event);
            }
        }
    }

    // Supervisor states
    StatesTableSTL c;
    StatesTableSTL rmtable;

    // f is a stack of states accessed in a dfs
    QStack<cldes_size_t> f;

    // Initialize f and ftable with the initial state
    f.push(virtualsys.init_state_);

    while (!f.isEmpty()) {
        auto const q = f.pop();

        if (!rmtable.contains(q) && !c.contains(q)) {
            // q = (qx, qy)
            auto const qx = q % aP.states_number_;
            auto const qy = q / aP.states_number_;

            auto const q_events =
              (aP.states_events_[qx] & aE.states_events_[qy]) |
              (aP.states_events_[qx] & virtualsys.only_in_0_) |
              (aE.states_events_[qy] & virtualsys.only_in_1_);

            auto const in_ncqx = p_non_contr_bit & aP.states_events_[qx];
            auto const in_ncqx_and_q = in_ncqx & q_events;

            if (in_ncqx_and_q != in_ncqx) {
                auto event = 0ul;
                auto event_it = in_ncqx_and_q ^ in_ncqx;
                while (event_it.any()) {
                    if (event_it.test(0)) {
                        RemoveBadStates(virtualsys,
                                        aP,
                                        aE,
                                        p_invgraph,
                                        e_invgraph,
                                        c,
                                        q,
                                        non_contr_bit,
                                        rmtable);
                    }
                    ++event;
                    event_it >>= 1ul;
                }
            } else {
                c.insert(q);

                virtualsys.transtriplet_.push_back(std::make_pair(
                  q, std::vector<std::pair<cldes_size_t, EventsBitArray>>()));

                auto event = 0ul;
                auto event_it = q_events;
                while (event_it.any()) {
                    if (event_it.test(0)) {
                        auto const fsqe = TransitionVirtual(aP, aE, q, event);

                        if (!rmtable.contains(fsqe)) {
                            if (!c.contains(fsqe)) {
                                f.push(fsqe);
                            }
                            virtualsys.transtriplet_.back().second.push_back(
                              std::make_pair(fsqe,
                                             EventsBitArray{ 1ul << event }));
                        }
                    }
                    ++event;
                    event_it >>= 1ul;
                }
            }
        }
    }

    rmtable.clear();

    // Swap new system states and sort it
    virtualsys.virtual_states_ = c.toList();
    qSort(virtualsys.virtual_states_);
    c.clear();

    // Make virtualsys a real sys
    SynchronizeStage2(virtualsys, aP, aE);

    virtualsys.Trim();

    return virtualsys;
}
