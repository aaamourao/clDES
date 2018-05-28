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
#include "des/transition_proxy.hpp"

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
        static_cast<cl_uint>(aSys0.device_graph_->size1()), asys0_private,
        aSys1.device_graph_->handle1().opencl_handle(),
        aSys1.device_graph_->handle2().opencl_handle(),
        aSys1.device_graph_->handle().opencl_handle(), asys1_private,
        result_dev.handle().opencl_handle(),
        static_cast<cl_uint>(result_dev.internal_size1())));

    // Copy device graph to host memory
    DESystemCL sync_sys(table_size, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));

    return sync_sys;
}
*/

DESystem op::Synchronize(DESystem &aSys0, DESystem &aSys1) {
    auto stage1 = SynchronizeStage1(aSys0, aSys1);

    return SynchronizeStage2(stage1, aSys0, aSys1);
}

op::StatesTableSTL op::SynchronizeStage1(DESystem const &aSys0,
                                         DESystem const &aSys1) {
    // Get the result and saves it on host memory
    StatesTableSTL states_table;

    for (auto ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
        for (auto ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
            states_table[ix1 * aSys0.states_number_ + ix0] =
                std::make_pair(ix0, ix1);
        }
    }

    return states_table;
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
        static_cast<cl_uint>(aSys0.device_graph_->size1()), asys0_private,
        aSys1.device_graph_->handle1().opencl_handle(),
        aSys1.device_graph_->handle2().opencl_handle(),
        aSys1.device_graph_->handle().opencl_handle(), asys1_private,
        result_dev.handle().opencl_handle(),
        static_cast<cl_uint>(result_dev.internal_size1())));

    // Copy device graph to host memory
    DESystemCL sync_sys(aTable->tsize, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));

    return sync_sys;
}
*/

DESystem op::SynchronizeStage2(op::StatesTableSTL const &aTable,
                               DESystem &aSys0, DESystem &aSys1) {
    std::set<cldes_size_t> markedstates_sync;

    DESystem result{aTable.size(), 0, markedstates_sync};

    auto sync_events = aSys0.events_ | aSys1.events_;

    for (auto q_ref = aTable.begin(); q_ref != aTable.end(); ++q_ref) {
        auto event = 0u;
        auto q = q_ref->second;
        auto sync_events_iter = sync_events;
        while (sync_events_iter != 0) {
            if (sync_events_iter[0]) {
                bool const is_in_p = aSys0.events_[event];
                bool const is_in_e = aSys1.events_[event];

                bool const is_in_x = (aSys0.states_events_[q.first])[event];
                bool const is_in_y = (aSys1.states_events_[q.second])[event];

                if ((is_in_x && is_in_y) || (is_in_x && !is_in_e) ||
                    (is_in_y && !is_in_p)) {
                    int xid = -1;
                    int yid;

                    auto index0 = std::distance(aTable.begin(), q_ref);

                    if (is_in_x) {
                        /*
                         * TODO: ublas::compressed_matrix<>::iterator1 does not
                         * have operators + and += implemented. Report bug and
                         * remove the following workaround.
                         */
                        auto p_row = aSys0.graph_.begin1();
                        for (auto i = 0; i < q.first; ++i) {
                            ++p_row;
                        }
                        // TODO: End of the workaround

                        for (auto elem = p_row.begin(); elem != p_row.end();
                             ++elem) {
                            if ((*elem)[event]) {
                                xid = elem.index2();
                                if (!is_in_e) {
                                    yid = q.second;
                                    auto index1_iter = aTable.find(
                                        yid * aSys0.states_number_ + xid);
                                    if (index1_iter != aTable.end()) {
                                        size_t index1 = std::distance(
                                            aTable.begin(), index1_iter);
                                        result(index0, index1) = event;
                                    }
                                    break;
                                }
                            }
                        }
                    }

                    if (is_in_y) {
                        /*
                         * TODO: ublas::compressed_matrix<>::iterator1 does not
                         * have operators + and += implemented. Report bug and
                         * remove the following workaround.
                         */
                        auto e_row = aSys1.graph_.begin1();
                        for (auto i = 0; i < q.second; ++i) {
                            ++e_row;
                        }
                        // TODO: End of the workaround

                        for (auto elem = e_row.begin(); elem != e_row.end();
                             ++elem) {
                            if ((*elem)[event]) {
                                if (!is_in_p) {
                                    xid = q.first;
                                }
                                yid = elem.index2();
                                auto index1_iter = aTable.find(
                                    yid * aSys0.states_number_ + xid);
                                if (index1_iter != aTable.end()) {
                                    size_t index1 = std::distance(
                                        aTable.begin(), index1_iter);
                                    result(index0, index1) = event;
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

    for (auto s0 : aSys0.marked_states_) {
        for (auto s1 : aSys1.marked_states_) {
            auto s_it = aTable.find(s1 * aSys0.states_number_ + s0);
            if (s_it != aTable.end()) {
                auto marked_state = std::distance(aTable.begin(), s_it);
                markedstates_sync.emplace(marked_state);
            }
        }
    }
    result.marked_states_ = markedstates_sync;

    auto init_state_sync_iter = aTable.find(
        aSys1.init_state_ * aSys0.states_number_ + aSys0.init_state_);
    auto init_state_sync = std::distance(aTable.begin(), init_state_sync_iter);
    result.init_state_ = init_state_sync;

    return result;
}

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

op::StatesTupleSTL op::TransitionVirtual(DESystem const &aSys0,
                                         DESystem const &aSys1,
                                         op::StatesTupleSTL const q,
                                         ScalarType const event) {
    bool const is_in_p = aSys0.events_[event];
    bool const is_in_e = aSys1.events_[event];

    bool const is_in_x = (aSys0.states_events_[q.first])[event];
    bool const is_in_y = (aSys1.states_events_[q.second])[event];

    int xid = -1;
    int yid;

    StatesTupleSTL ret;

    if (is_in_x) {
        /*
         * TODO: ublas::compressed_matrix<>::iterator1 does not
         * have operators + and += implemented. Report bug and
         * remove the following workaround.
         */
        auto p_row = aSys0.graph_.begin1();
        for (auto i = 0; i < q.first; ++i) {
            ++p_row;
        }
        // TODO: End of the workaround

        for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
            if ((*elem)[event]) {
                xid = elem.index2();
                if (!is_in_e) {
                    yid = q.second;
                    ret.first = xid;
                    ret.second = yid;
                    return ret;
                }
            }
        }
    }

    if (is_in_y) {
        /*
         * TODO: ublas::compressed_matrix<>::iterator1 does not have
         * operators + and += implemented. Report bug and remove the
         * following workaround.
         */
        auto e_row = aSys1.graph_.begin1();
        for (auto i = 0; i < q.second; ++i) {
            ++e_row;
        }
        // TODO: End of the workaround

        for (auto elem = e_row.begin(); elem != e_row.end(); ++elem) {
            if ((*elem)[event]) {
                if (!is_in_p) {
                    xid = q.first;
                }
                yid = elem.index2();
                ret.first = xid;
                ret.second = yid;
                return ret;
            }
        }
    }

    return ret;
}

bool op::ExistTransitionReal(DESystem const &aSys, cldes_size_t const &x,
                             ScalarType const &event) {
    bool const is_in_x = (aSys.states_events_[x])[event];
    return is_in_x;
}

template <class EventsType, class GraphType>
static op::StatesTableSTL __TransitionVirtualInv(EventsType const &aEventsP,
                                                 EventsType const &aEventsE,
                                                 GraphType const &aInvGraphP,
                                                 GraphType const &aInvGraphE,
                                                 op::StatesTupleSTL const q,
                                                 ScalarType const event) {
    bool const is_in_p = aEventsP[event];
    bool const is_in_e = aEventsE[event];

    op::StatesTableSTL ret;
    if (!is_in_p && !is_in_e) {
        return ret;
    }

    /*
     * TODO: ublas::compressed_matrix<>::iterator1 does not have
     * operators + and += implemented. Report bug and remove the
     * following workaround.
     */
    auto p_row = aInvGraphP.begin1();
    for (auto i = 0; i < q.first; ++i) {
        ++p_row;
    }
    auto e_row = aInvGraphE.begin1();
    for (auto i = 0; i < q.second; ++i) {
        ++e_row;
    }
    // TODO: End of the workaround

    if (is_in_p) {
        for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
            if ((*elem)[event]) {
                if (!is_in_e) {
                    auto key = q.second * aInvGraphP.size1() + elem.index2();
                    ret[key] = std::make_pair(elem.index2(), q.second);
                } else {
                    for (auto elemg2 = e_row.begin(); elemg2 != e_row.end();
                         ++elemg2) {
                        if ((*elemg2)[event]) {
                            auto key = elemg2.index2() * aInvGraphP.size1() +
                                       elem.index2();
                            ret[key] =
                                std::make_pair(elem.index2(), elemg2.index2());
                        }
                    }
                }
            }
        }
    } else {  // Is only on e: !is_in_p && is_in_e
        for (auto elemg2 = e_row.begin(); elemg2 != e_row.end(); ++elemg2) {
            if ((*elemg2)[event]) {
                auto key = elemg2.index2() * aInvGraphP.size1() + q.first;
                ret[key] = std::make_pair(q.first, elemg2.index2());
            }
        }
    }

    return ret;
}

template <class EventsType, class GraphType>
static void __RemoveBadStates(EventsType const &aEventsP,
                              EventsType const &aEventsE,
                              GraphType const &aInvGraphP,
                              GraphType const &aInvGraphE,
                              op::StatesTableSTL &C,
                              op::StatesTupleSTL const &q,
                              std::unordered_set<ScalarType> const &s_non_contr,
                              op::StatesTableSTL &aRemovedStates) {
    op::StatesTableSTL f;
    f[q.second * aInvGraphP.size1() + q.first] = q;
    aRemovedStates[q.second * aInvGraphP.size1() + q.first] = q;

    while (f.size() != 0) {
        op::StatesTupleSTL x = std::get<1>(*(f.begin()));
        C.erase(x.second * aInvGraphP.size1() + x.first);
        f.erase(f.begin());
        for (auto e : s_non_contr) {
            auto finv = __TransitionVirtualInv(aEventsP, aEventsE, aInvGraphP,
                                               aInvGraphE, x, e);
            if (finv.size() != 0) {
                f.insert(finv.begin(), finv.end());
                aRemovedStates.insert(finv.begin(), finv.end());
            }
        }
    }

    return;
}

DESystem op::SupervisorSynth(DESystem &aP, DESystem &aE,
                             std::unordered_set<ScalarType> const &non_contr) {
    auto s_events = aP.events_ | aE.events_;

    auto const p_inversedgraph = boost::numeric::ublas::trans(aP.graph_);
    auto const e_inversedgraph = boost::numeric::ublas::trans(aE.graph_);

    auto states_table = SynchronizeStage1(aP, aE);

    auto init_state_sync_key =
        aE.init_state_ * aP.states_number_ + aP.init_state_;
    auto init_state_sync_iter = states_table.find(init_state_sync_key);

    StatesTableSTL f;

    if (init_state_sync_iter != states_table.end()) {
        f[init_state_sync_key] = std::make_pair(aP.init_state_, aE.init_state_);
    }

    StatesTableSTL s_states;
    StatesTableSTL removed_states;

    while (f.size() != 0) {
        auto q = std::get<1>(*(f.begin()));
        s_states[q.second * aP.states_number_ + q.first] = q;
        f.erase(f.begin());
        auto event = 0u;
        auto s_events_iter = s_events;
        while (s_events_iter != 0) {
            if (s_events_iter[0]) {
                bool const is_non_contr =
                    non_contr.find(event) != non_contr.end();
                auto is_fp = ExistTransitionReal(aP, q.first, event);
                auto is_fs_qevent = ExistTransitionVirtual(aP, aE, q, event);
                StatesTupleSTL fs_qevent;
                if (is_fs_qevent) {
                    fs_qevent = TransitionVirtual(aP, aE, q, event);
                    if (removed_states.find(
                            fs_qevent.second * aP.states_number_ +
                            fs_qevent.first) != removed_states.end()) {
                        is_fs_qevent = false;
                    }
                }
                if (is_non_contr && !is_fs_qevent && is_fp) {
                    __RemoveBadStates(aP.events_, aE.events_, p_inversedgraph,
                                      e_inversedgraph, s_states, q, non_contr,
                                      removed_states);
                    break;
                } else if (is_fs_qevent) {
                    auto key =
                        fs_qevent.second * aP.states_number_ + fs_qevent.first;
                    auto is_in_f = f.find(key) != f.end();
                    auto is_in_s_states = s_states.find(key) != s_states.end();
                    if (!is_in_s_states && !is_in_f) {
                        f[key] = fs_qevent;
                    }
                }
            }
            ++event;
            s_events_iter >>= 1u;
        }
    }

    auto supervisor = SynchronizeStage2(s_states, aP, aE);

    for (auto st : states_table) {
        auto s = st.second;
        std::cout << "(" << s.first << ", " << s.second << "): ";
        for (auto e = 0; e < 64; ++e) {
            auto finv = __TransitionVirtualInv(
                aP.events_, aE.events_, p_inversedgraph, e_inversedgraph, s, e);
            for (auto qt : finv) {
                auto q = qt.second;
                std::cout << "((" << q.first << ", " << q.second << "), " << e
                          << ") ";
            }
        }
        std::cout << std::endl;
    }

    return supervisor.Trim();
}
