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

#include "operations/operations.hpp"
#include <algorithm>
#include <cmath>
// #include "backend/oclbackend.hpp"
#include "des/desystem.hpp"
// #include "des/desystemcl.hpp"
#include <eigen3/Eigen/Sparse>
#include "des/transition_proxy.hpp"

#include <iostream>

using namespace cldes;

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

DESystem op::Synchronize(DESystem &aSys0, DESystem &aSys1) {
    auto syncsys = SynchronizeStage1(aSys0, aSys1);
    SynchronizeStage2(syncsys, aSys0, aSys1);
    return syncsys;
}

DESystem op::SynchronizeStage1(DESystem const &aSys0, DESystem const &aSys1) {
    auto in_both = aSys0.events_ & aSys1.events_;
    auto only_in_0 = aSys0.events_ ^ in_both;
    auto only_in_1 = aSys1.events_ ^ in_both;

    // New system params
    auto estimated_transitions = 0;
    DESystem::EventsSet events = 0ull;
    DESystem::StatesEventsTable states_events;
    DESystem::StatesEventsTable inv_states_events;
    states_events.reserve(aSys0.states_number_ * aSys1.states_number_);
    inv_states_events.reserve(aSys0.states_number_ * aSys1.states_number_);

    // Calculate params
    for (auto ix0 = 0ul; ix0 < aSys0.states_number_; ++ix0) {
        for (auto ix1 = 0ul; ix1 < aSys1.states_number_; ++ix1) {
            auto key = ix1 * aSys0.states_number_ + ix0;

            states_events[key] =
                (aSys0.states_events_[ix0] & aSys1.states_events_[ix1]) |
                (aSys0.states_events_[ix0] & only_in_0) |
                (aSys1.states_events_[ix1] & only_in_1);
            inv_states_events[key] =
                (aSys0.inv_states_events_[ix0] &
                 aSys1.inv_states_events_[ix1]) |
                (aSys0.inv_states_events_[ix0] & only_in_0) |
                (aSys1.inv_states_events_[ix1] & only_in_1);

            estimated_transitions += states_events[key].count();
            events |= states_events[key];
        }
    }

    // Calculate new marked states
    DESystem::StatesSet marked_states;
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            marked_states.insert(s1 * aSys0.states_number_ + s0);
        }
    }

    // Create new system without transitions
    DESystem virtualsys{
        aSys0.states_number_ * aSys1.states_number_,
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_,
        marked_states};

    // Set private params
    virtualsys.events_ = events;
    virtualsys.states_events_ = states_events;
    virtualsys.inv_states_events_ = inv_states_events;

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

using RowIterator = Eigen::InnerIterator<DESystem::GraphHostData const>;

void op::SynchronizeStage2(DESystem &aVirtualSys, DESystem const &aSys0,
                           DESystem const &aSys1) {
    /*
    std::vector<std::pair(StatesTupleSTL, cldes_size_t)> stable;
    for (auto sit = aVirtualSys.states_events_.begin();
         sit != aVirtualSys.states_events_.end(); ++sit) {
        cldes_size_t pos =
            std::distance(aVirtualStates.states_events_.begin(), sit);
        stable.push_back(
            std::make_pair(std::make_pair(sit->first % aSys0.states_number_,
                                          sit->first / aSys0.states_number_),
                           pos));
    }
    */

    // Make it const to avoid eventual re-hashes due to the insertion of
    // transitions
    Eigen::VectorXi transitions_number(aVirtualSys.states_events_.size());
    std::map<cldes_size_t, EventsBitArray> virtualse;
    for (auto st = aVirtualSys.states_events_.constBegin();
         st != aVirtualSys.states_events_.constEnd(); ++st) {
        transitions_number(virtualse.size()) = st.value().count();
        virtualse[st.key()] = st.value();
    }

    // Delete states_events_, it will be remapped
    aVirtualSys.states_events_.clear();
    aVirtualSys.inv_states_events_.clear();

    // Resize adj matrices if necessary
    if (aVirtualSys.states_number_ != virtualse.size()) {
        aVirtualSys.graph_.resize(virtualse.size(), virtualse.size());
        aVirtualSys.bit_graph_.resize(virtualse.size(), virtualse.size());
        aVirtualSys.states_number_ = virtualse.size();

        // Initialize bit graph with Identity
        aVirtualSys.bit_graph_.setIdentity();
    }

    // Reserve space for transitions
    aVirtualSys.graph_.reserve(transitions_number);
    aVirtualSys.bit_graph_.reserve(transitions_number);

    // Calculate transitions
    for (auto q_ref = virtualse.begin(); q_ref != virtualse.end(); ++q_ref) {
        auto q = std::make_pair(q_ref->first % aSys0.states_number_,
                                q_ref->first / aSys0.states_number_);
        auto sync_events_iter = aVirtualSys.events_;
        auto event = 0u;
        while (sync_events_iter != 0) {
            if (sync_events_iter[0]) {
                if (q_ref->second[event]) {
                    auto index0 = std::distance(virtualse.begin(), q_ref);

                    int xto;
                    int yto;

                    bool const is_in_p = aSys0.events_[event];
                    bool const is_in_e = aSys1.events_[event];

                    if (is_in_p) {
                        for (RowIterator pe(aSys0.graph_, q.first); pe; ++pe) {
                            if (pe.value()[event]) {
                                xto = pe.col();
                                if (!is_in_e) {
                                    yto = q.second;
                                    auto index1_iter = virtualse.find(
                                        yto * aSys0.states_number_ + xto);
                                    if (index1_iter != virtualse.end()) {
                                        size_t index1 = std::distance(
                                            virtualse.begin(), index1_iter);
                                        aVirtualSys(index0, index1) = event;
                                    }
                                    break;
                                }
                            }
                        }
                    }

                    if (is_in_e) {
                        for (RowIterator ee(aSys1.graph_, q.second); ee; ++ee) {
                            if (ee.value()[event]) {
                                if (!is_in_p) {
                                    xto = q.first;
                                }
                                yto = ee.col();
                                auto index1_iter = virtualse.find(
                                    yto * aSys0.states_number_ + xto);
                                if (index1_iter != virtualse.end()) {
                                    size_t index1 = std::distance(
                                        virtualse.begin(), index1_iter);
                                    aVirtualSys(index0, index1) = event;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            ++event;
            sync_events_iter >>= 1u;
        }
    }

    // Remove aditional space
    aVirtualSys.graph_.makeCompressed();
    aVirtualSys.bit_graph_.makeCompressed();

    // Remap marked states
    aVirtualSys.marked_states_.clear();
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            auto marked_state_it =
                virtualse.find(s1 * aSys0.states_number_ + s0);
            size_t marked_state =
                std::distance(virtualse.begin(), marked_state_it);
            aVirtualSys.marked_states_.insert(marked_state);
        }
    }

    // Remap initial state
    auto init_state_sync_iter = virtualse.find(
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_);
    size_t init_state = std::distance(virtualse.begin(), init_state_sync_iter);
    aVirtualSys.init_state_ = init_state;
}

/*
bool op::ExistTransitionVirtual(DESystem const &aSys0, DESystem const &aSys1,
                                op::StatesTupleSTL const q,
                                ScalarType const event) {
    bool const is_in_p = aSys0.events_[event];
    bool const is_in_e = aSys1.events_[event];

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

op::StatesTupleSTL op::TransitionVirtual(DESystem const &aSys0,
                                         DESystem const &aSys1,
                                         cldes_size_t const &q,
                                         ScalarType const &event) {
    bool const is_in_p = aSys0.events_[event];
    bool const is_in_e = aSys1.events_[event];

    auto qx = q % aSys0.states_number_;
    auto qy = q / aSys0.states_number_;

    int xid;
    int yid;

    StatesTupleSTL ret;

    if (is_in_p) {
        for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
            if (pe.value()[event]) {
                xid = pe.col();
                if (!is_in_e) {
                    yid = qy;
                    ret.first = xid;
                    ret.second = yid;
                    return ret;
                }
                break;
            }
        }
    }

    if (is_in_e) {
        for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
            if (ee.value()[event]) {
                if (!is_in_p) {
                    xid = qx;
                }
                yid = ee.col();
                ret.first = xid;
                ret.second = yid;
                return ret;
            }
        }
    }
    return ret;
}

// This function assumes that there is an inverse transition.
template <class EventsType>
static op::StatesTableSTL __TransitionVirtualInv(
    EventsType const &aEventsP, EventsType const &aEventsE,
    op::GraphType const &aInvGraphP, op::GraphType const &aInvGraphE,
    cldes_size_t const &q, ScalarType const &event) {
    auto qx = q % aInvGraphP.rows();
    auto qy = q / aInvGraphP.rows();

    bool const is_in_p = aEventsP[event];
    bool const is_in_e = aEventsE[event];

    op::StatesTableSTL ret;
    ret.reserve(aEventsP.size());

    if (is_in_p && is_in_e) {
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value()[event]) {
                for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
                    if (ee.value()[event]) {
                        auto key = ee.col() * aInvGraphP.rows() + pe.col();
                        ret.insert(key);
                    }
                }
            }
        }
    } else if (is_in_p) {  // Is only in p: is_in_p && !is_in_e
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value()[event]) {
                ret.insert(qy * aInvGraphP.rows() + pe.col());
            }
        }
    } else {  // Is only in e: !is_in_p && is_in_e
        for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
            if (ee.value()[event]) {
                ret.insert(ee.col() * aInvGraphP.rows() + qx);
            }
        }
    }

    return ret;
}

void op::RemoveBadStates(DESystem &aVirtualSys, DESystem const &aP,
                         DESystem const &aE, op::GraphType const &aInvGraphP,
                         op::GraphType const &aInvGraphE, op::StatesTableSTL &C,
                         op::StatesTableSTL &fs, cldes_size_t const &q,
                         QSet<ScalarType> const &s_non_contr) {
    StatesTableSTL f;
    f.reserve(aVirtualSys.states_number_);
    f.insert(q);

    while (f.size() != 0) {
        cldes_size_t const x = *(f.begin());

        C.remove(x);
        f.remove(x);
        fs.remove(x);

        auto events_iter = aVirtualSys.inv_states_events_[x];
        auto event = 0u;
        while (events_iter != 0) {
            if (events_iter[0]) {
                auto finv = __TransitionVirtualInv(
                    aP.events_, aE.events_, aInvGraphP, aInvGraphE, x, event);
                auto is_non_contr =
                    s_non_contr.find(event) != s_non_contr.end();

                if (is_non_contr) {
                    for (auto s : finv) {
                        if (aVirtualSys.states_events_.find(s) !=
                            aVirtualSys.states_events_.end()) {
                            // Add transitions to remove set
                            f.insert(s);
                        }
                    }
                } else {
                    // Remove transitions to removed state
                    for (auto s : finv) {
                        if (aVirtualSys.states_events_.find(s) !=
                            aVirtualSys.states_events_.end()) {
                            aVirtualSys.states_events_[s][event] = false;
                        }
                    }
                }
            }
            ++event;
            events_iter >>= 1u;
        }
        // Remove states from hash table
        aVirtualSys.states_events_.remove(x);
        aVirtualSys.inv_states_events_.remove(x);
    }

    return;
}

DESystem op::SupervisorSynth(DESystem const &aP, DESystem const &aE,
                             QSet<ScalarType> const &non_contr) {
    DESystem::GraphHostData const p_invgraph = aP.graph_.transpose();
    DESystem::GraphHostData const e_invgraph = aE.graph_.transpose();

    auto virtualsys = SynchronizeStage1(aP, aE);

    QSet<ScalarType> s_non_contr;
    for (auto eu : non_contr) {
        if (virtualsys.events_[eu]) {
            s_non_contr.insert(eu);
        }
    }

    StatesTableSTL c;
    StatesTableSTL f;
    c.reserve(virtualsys.states_number_);
    f.reserve(virtualsys.states_number_);

    f.insert(virtualsys.init_state_);

    while (f.size() != 0) {
        auto q = *(f.begin());
        f.remove(q);
        c.insert(q);

        // q = (qx, qy)
        auto qx = q % aP.states_number_;

        auto event = 0u;
        auto s_events_iter = virtualsys.events_;
        while (s_events_iter != 0) {
            if (s_events_iter[0]) {
                bool const is_non_contr =
                    s_non_contr.find(event) != s_non_contr.end();
                auto is_there_fp = aP.states_events_[qx][event];
                auto is_there_fsqe = virtualsys.states_events_[q][event];

                if (is_non_contr && !is_there_fsqe && is_there_fp) {
                    RemoveBadStates(virtualsys, aP, aE, p_invgraph, e_invgraph,
                                    c, f, q, s_non_contr);
                    break;
                } else if (is_there_fsqe) {
                    auto fs_qevent = TransitionVirtual(aP, aE, q, event);
                    auto fsqe_key =
                        fs_qevent.second * aP.states_number_ + fs_qevent.first;

                    auto is_in_f = f.find(fsqe_key) != f.end();
                    auto is_in_c = c.find(fsqe_key) != c.end();

                    if (!is_in_c && !is_in_f) {
                        f.insert(fsqe_key);
                    }
                }
            }
            ++event;
            s_events_iter >>= 1u;
        }
    }

    auto const syscopy = virtualsys.states_events_;
    for (auto st = syscopy.constBegin(); st != syscopy.constEnd(); ++st) {
        if (c.find(st.key()) == c.end()) {
            virtualsys.states_events_.remove(st.key());
        }
    }

    // Make virtualsys a real sys
    SynchronizeStage2(virtualsys, aP, aE);

    /*
    DESystem::StatesSet bla;
    DESystem test{6, 0, bla};
    test(0, 1) = 1;
    test(0, 2) = 0;
    test(0, 2) = 2;
    test(1, 1) = 1;
    test(1, 2) = 0;
    test(1, 3) = 2;
    test(2, 3) = 1;
    test(2, 4) = 0;
    test(3, 3) = 1;
    test(3, 4) = 0;
    test(4, 0) = 2;
    test(4, 4) = 0;
    test(5, 1) = 2;
    test(5, 4) = 0;

    cldes::DESystem g1{3, 0, bla};

    g1(0, 0) = 0;
    g1(0, 2) = 2;
    g1(1, 0) = 0;
    g1(1, 1) = 1;
    g1(2, 1) = 0;
    g1(2, 1) = 2;
    g1(2, 2) = 1;

    cldes::DESystem g2{2, 0, bla};

    g2(0, 0) = 1;
    g2(0, 1) = 0;
    g2(1, 0) = 1;
    g2(1, 1) = 0;

    // for (auto st : states_table) {
    //    auto s = st.second;
    // std::cout << "(" << s.first << ", " << s.second << "): ";
    for (auto s0 = g1.states_events_.begin(); s0 != g1.states_events_.end();
         ++s0) {
        for (auto s1 = g2.states_events_.begin(); s1 != g2.states_events_.end();
             ++s1) {
            auto i = s1->first * g1.states_number_ + s0->first;
            std::cout << i << ": ";
            // for (auto e = 0; e < 64; ++e) {
            auto finv = __TransitionVirtualInv(g1.events_, g2.events_,
                                               ublas::trans(g1.graph_),
                                               ublas::trans(g2.graph_), i, 0);
            for (auto qt : finv) {
                // auto qx = qt % 3;
                // auto qy = qt / 3;
                // std::cout << "((" << qx << ", " << qy << "), 0) ";
                std::cout << qt << ":0 ";
            }
            //    }
            std::cout << std::endl;
            //}
        }
    }
*/
    virtualsys.Trim();

    return virtualsys;
}
