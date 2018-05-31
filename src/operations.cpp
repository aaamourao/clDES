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
#include <qt5/QtCore/QVector>
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
    auto const in_both = aSys0.events_ & aSys1.events_;
    auto const only_in_0 = aSys0.events_ ^ in_both;
    auto const only_in_1 = aSys1.events_ ^ in_both;

    // New system params
    DESystem::StatesEventsTable states_events;
    DESystem::StatesEventsTable inv_states_events;
    states_events.reserve(aSys0.states_number_ * aSys1.states_number_);
    inv_states_events.reserve(aSys0.states_number_ * aSys1.states_number_);

    // Calculate params
    for (auto ix0 = 0ul; ix0 < aSys0.states_number_; ++ix0) {
        for (auto ix1 = 0ul; ix1 < aSys1.states_number_; ++ix1) {
            auto const key = ix1 * aSys0.states_number_ + ix0;

            states_events[key] = (aSys0.states_events_.value(ix0) &
                                  aSys1.states_events_.value(ix1)) |
                                 (aSys0.states_events_.value(ix0) & only_in_0) |
                                 (aSys1.states_events_.value(ix1) & only_in_1);
            inv_states_events[key] =
                (aSys0.inv_states_events_.value(ix0) &
                 aSys1.inv_states_events_.value(ix1)) |
                (aSys0.inv_states_events_.value(ix0) & only_in_0) |
                (aSys1.inv_states_events_.value(ix1) & only_in_1);
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
    virtualsys.events_ = aSys0.events_ | aSys1.events_;
    states_events.swap(virtualsys.states_events_);
    inv_states_events.swap(virtualsys.inv_states_events_);

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

    // Calculate sparcity pattern and create efficient data structures for
    // searching states
    QList<cldes_size_t> states;
    QHash<cldes_size_t, cldes_size_t> statesmap;
    Eigen::VectorXi sparcitypattern(aVirtualSys.states_events_.size());

    aVirtualSys.states_events_.keys().swap(states);
    qSort(states);

    statesmap.reserve(states.size());

    auto pos = 0ul;
    foreach (cldes_size_t s, states) {
        sparcitypattern.coeffRef(pos) =
            aVirtualSys.states_events_.value(s).count();
        statesmap[s] = pos;
        ++pos;
    }

    // Copy and delete states_events_, it will be remapped
    auto const virtualse = aVirtualSys.states_events_;
    aVirtualSys.states_events_.clear();
    aVirtualSys.inv_states_events_.clear();

    // Resize adj matrices if necessary
    if (static_cast<long>(aVirtualSys.states_number_) != states.size()) {
        aVirtualSys.graph_.resize(states.size(), states.size());
        aVirtualSys.bit_graph_.resize(states.size(), states.size());
        aVirtualSys.states_number_ = states.size();
    }

    // Reserve space for transitions
    aVirtualSys.graph_.reserve(sparcitypattern);
    aVirtualSys.bit_graph_.reserve(sparcitypattern);
    aVirtualSys.states_events_.reserve(states.size());
    aVirtualSys.inv_states_events_.reserve(states.size());

    // Initialize bit graph with Identity
    aVirtualSys.bit_graph_.setIdentity();

    // Calculate transitions
    foreach (cldes_size_t s, states) {
        auto const qx = s % aSys0.states_number_;
        auto const qy = s / aSys0.states_number_;

        auto q_events = virtualse.value(s);

        auto event = 0ul;
        while (q_events.any()) {
            if (q_events.test(0)) {
                int xto;
                int yto;

                bool const is_in_p = aSys0.events_.test(event);
                bool const is_in_e = aSys1.events_.test(event);

                if (is_in_p) {
                    for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
                        if (pe.value().test(event)) {
                            xto = pe.col();
                            if (!is_in_e) {
                                yto = qy;
                                break;
                            }
                        }
                    }
                }

                if (is_in_e) {
                    for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
                        if (ee.value().test(event)) {
                            if (!is_in_p) {
                                xto = qx;
                            }
                            yto = ee.col();
                            break;
                        }
                    }
                }

                auto const key = yto * aSys0.states_number_ + xto;
                if (statesmap.contains(key)) {
                    aVirtualSys(statesmap.value(s), statesmap.value(key)) =
                        event;
                }
            }
            ++event;
            q_events >>= 1ul;
        }
    }

    // Remove aditional space
    aVirtualSys.graph_.makeCompressed();
    aVirtualSys.bit_graph_.makeCompressed();

    // Remap marked states
    aVirtualSys.marked_states_.clear();
    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            auto const key = s1 * aSys1.states_number_ + s0;
            if (statesmap.contains(key)) {
                aVirtualSys.marked_states_.insert(
                    statesmap.value(s1 * aSys1.states_number_ + s0));
            }
        }
    }

    // It only works for init_state = 0;
    aVirtualSys.init_state_ =
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_;
    /*
    // Remap initial state
    auto const init_state_key =
        aSys1.init_state_ + aSys0.states_number_ + aSys0.init_state_;
    if (statesmap.find(init_state_key) != statesmap.end()) {
        aVirtualSys.init_state_ = statesmap[init_state_key];
    } else {
        // Set init state to an invalid state
        aVirtualSys.init_state_ = aVirtualSys.states_number_;
    }
    */
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
    bool const is_in_p = aSys0.events_.test(event);
    bool const is_in_e = aSys1.events_.test(event);

    auto const qx = q % aSys0.states_number_;
    auto const qy = q / aSys0.states_number_;

    int xid;
    int yid;

    if (is_in_p) {
        for (RowIterator pe(aSys0.graph_, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                xid = pe.col();
                if (!is_in_e) {
                    yid = qy;
                }
                break;
            }
        }
    }

    if (is_in_e) {
        for (RowIterator ee(aSys1.graph_, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                if (!is_in_p) {
                    xid = qx;
                }
                yid = ee.col();
                break;
            }
        }
    }
    StatesTupleSTL ret;
    ret.first = xid;
    ret.second = yid;
    return ret;
}

using StatesArray = QVector<cldes_size_t>;

// This function assumes that there is an inverse transition.
template <class EventsType>
static StatesArray __TransitionVirtualInv(EventsType const &aEventsP,
                                          EventsType const &aEventsE,
                                          op::GraphType const &aInvGraphP,
                                          op::GraphType const &aInvGraphE,
                                          cldes_size_t const &q,
                                          ScalarType const &event) {
    auto const qx = q % aInvGraphP.rows();
    auto const qy = q / aInvGraphP.rows();

    bool const is_in_p = aEventsP.test(event);
    bool const is_in_e = aEventsE.test(event);

    StatesArray ret;
    // ret.reserve(cldes::g_max_events);

    auto const p_size = aInvGraphP.rows();

    if (is_in_p && is_in_e) {
        StatesArray pstates;
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                pstates.push_back(pe.col());
            }
        }
        for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                foreach (cldes_size_t sp, pstates) {
                    ret.push_back(ee.col() * p_size + sp);
                }
            }
        }
    } else if (is_in_p) {  // Is only in p: is_in_p && !is_in_e
        for (RowIterator pe(aInvGraphP, qx); pe; ++pe) {
            if (pe.value().test(event)) {
                ret.push_back(qy * p_size + pe.col());
            }
        }
    } else {  // Is only in e: !is_in_p && is_in_e
        for (RowIterator ee(aInvGraphE, qy); ee; ++ee) {
            if (ee.value().test(event)) {
                ret.push_back(ee.col() * p_size + qx);
            }
        }
    }

    return ret;
}

void op::RemoveBadStates(DESystem &aVirtualSys, DESystem const &aP,
                         DESystem const &aE, op::GraphType const &aInvGraphP,
                         op::GraphType const &aInvGraphE,
                         QHash<cldes_size_t, EventsBitArray> &C,
                         op::StatesStack &fs, cldes_size_t const &q,
                         QSet<ScalarType> const &s_non_contr) {
    StatesStack f;
    f.push(q);

    while (!f.isEmpty()) {
        cldes_size_t const x = f.pop();

        C.remove(x);
        fs.removeOne(x);

        auto const q_events = aVirtualSys.inv_states_events_.value(x);

        foreach (ScalarType event, s_non_contr) {
            if (q_events.test(event)) {
                auto const finv = __TransitionVirtualInv(
                    aP.events_, aE.events_, aInvGraphP, aInvGraphE, x, event);

                foreach (cldes_size_t s, finv) {
                    if (aVirtualSys.states_events_.contains(s)) {
                        f.push(s);
                    }
                }
            }
        }
        aVirtualSys.states_events_.remove(x);
    }

    return;
}

DESystem op::SupervisorSynth(DESystem const &aP, DESystem const &aE,
                             QSet<ScalarType> const &non_contr) {
    DESystem::GraphHostData const p_invgraph = aP.graph_.transpose();
    DESystem::GraphHostData const e_invgraph = aE.graph_.transpose();

    auto virtualsys = SynchronizeStage1(aP, aE);

    QSet<ScalarType> s_non_contr = non_contr;
    EventsBitArray non_contr_bit;
    foreach (ScalarType event, non_contr) {
        if (!virtualsys.events_.test(event)) {
            s_non_contr.remove(event);
        } else {
            non_contr_bit.set(event);
        }
    }

    DESystem::StatesEventsTable c;
    c.reserve(virtualsys.states_number_ * 3 / 100);

    StatesStack f;
    f.push(virtualsys.init_state_);

    while (!f.isEmpty()) {
        auto const q = f.pop();
        c[q] = virtualsys.states_events_.value(q);

        // q = (qx, qy)
        auto const qx = q % aP.states_number_;

        auto event = 0ul;
        auto event_it = virtualsys.states_events_.value(q) | non_contr_bit;
        while (event_it.any()) {
            if (event_it.test(0)) {
                auto const is_there_fsqe =
                    virtualsys.states_events_.value(q).test(event);

                if (s_non_contr.contains(event) && !is_there_fsqe &&
                    aP.states_events_.value(qx).test(event)) {
                    RemoveBadStates(virtualsys, aP, aE, p_invgraph, e_invgraph,
                                    c, f, q, s_non_contr);
                    break;
                } else if (is_there_fsqe) {
                    auto const fs_qevent = TransitionVirtual(aP, aE, q, event);
                    auto const fsqe_key =
                        fs_qevent.second * aP.states_number_ + fs_qevent.first;

                    if (!c.contains(fsqe_key) && !f.contains(fsqe_key)) {
                        f.push(fsqe_key);
                    }
                }
            }
            ++event;
            event_it >>= 1ul;
        }
    }

    c.swap(virtualsys.states_events_);
    virtualsys.states_events_.squeeze();

    // Make virtualsys a real sys
    SynchronizeStage2(virtualsys, aP, aE);

    virtualsys.Trim();

    return virtualsys;
}
