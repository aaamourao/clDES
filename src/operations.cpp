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
#include "backend/oclbackend.hpp"
#include "des/desystem.hpp"
#include "des/desystemcl.hpp"

using namespace cldes;

cldes_size_t op::TablePos_(cldes_size_t const &aG0Pos,
                           cldes_size_t const &aG1Pos,
                           cldes_size_t const &aG0NStates) {
    return aG1Pos * aG0NStates + aG0Pos;
}

float op::CalcGCD_(float aG0, float aG1) {
    if (aG0 == 0 || aG1 == 0) {
        return 1.0f;
    }

    float mod;

    while (aG1 != 0) {
        mod = std::fmod(aG0, aG1);
        aG0 = aG1;
        aG1 = mod;
    }

    return aG0;
}

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

DESystem op::Synchronize(DESystem &aSys0, DESystem &aSys1) {
    // Sync stage 1
    StatesTableSTL states_tuple;
    for (unsigned int ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
        for (unsigned int ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
            states_tuple.push_back(std::make_tuple(ix0, ix1));
        }
    }

    // Sync stage 2
    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    std::set<float> sync_events;
    std::set_union(aSys0.events_.begin(), aSys0.events_.end(),
                   aSys1.events_.begin(), aSys1.events_.end(),
                   std::inserter(sync_events, sync_events.begin()));

    DESystem::GraphHostData result(states_tuple.size(), states_tuple.size());

    for (auto q : states_tuple) {
        for (auto event : sync_events) {
            bool const is_in_p =
                std::find(aSys0.events_.begin(), aSys0.events_.end(), event) !=
                aSys0.events_.end();
            bool const is_in_e =
                std::find(aSys1.events_.begin(), aSys1.events_.end(), event) !=
                aSys1.events_.end();
            if (!is_in_p && !is_in_e) {
                break;
            }

            int xid = -1;
            int yid;

            if (is_in_p) {
                /*
                 * TODO: ublas::compressed_matrix<>::iterator1 does not have
                 * operators + and += implemented. Report bug and remove the
                 * following workaround.
                 */
                auto p_row = aSys0.graph_.begin1();
                for (auto i = 0; i < std::get<0>(q); ++i) {
                    ++p_row;
                }
                // TODO: End of the workaround

                for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
                    if (std::fmod(*elem, event) == 0.0f) {
                        xid = elem.index2();
                        if (!is_in_e) {
                            yid = std::get<1>(q);
                            auto index0 =
                                std::get<1>(q) * aSys0.states_number_ +
                                std::get<0>(q);
                            auto index1 = yid * aSys0.states_number_ + xid;
                            if (result(index0, index1) != 0) {
                                result(index0, index1) *= event;
                            } else {
                                result(index0, index1) = event;
                            }
                            break;
                        }
                    }
                }
            }

            /*
             * TODO: ublas::compressed_matrix<>::iterator1 does not have
             * operators + and += implemented. Report bug and remove the
             * following workaround.
             */
            auto e_row = aSys1.graph_.begin1();
            for (auto i = 0; i < std::get<1>(q); ++i) {
                ++e_row;
            }
            // TODO: End of the workaround

            if ((is_in_e && !is_in_p) || (is_in_e && (xid != -1))) {
                for (auto elem = e_row.begin(); elem != e_row.end(); ++elem) {
                    if (std::fmod(*elem, event) == 0.0f) {
                        if (!is_in_p) {
                            xid = std::get<0>(q);
                        }
                        yid = elem.index2();
                        auto index0 = std::get<1>(q) * aSys0.states_number_ +
                                      std::get<0>(q);
                        auto index1 = yid * aSys0.states_number_ + xid;
                        if (result(index0, index1) != 0) {
                            result(index0, index1) *= event;
                        } else {
                            result(index0, index1) = event;
                        }
                        break;
                    }
                }
            }
        }
    }

    // Copy device graph to host memory
    DESystem sync_sys(result, states_tuple.size(), initstate_sync,
                      markedstates_sync);
    sync_sys.events_ = sync_events;

    return sync_sys;
}

op::StatesTableSTL op::SynchronizeStage1(DESystem const &aSys0,
                                         DESystem const &aSys1) {
    // Get the result and saves it on host memory
    StatesTableSTL states_table;

    for (auto ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
        for (auto ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
            states_table.push_back(std::make_tuple(ix0, ix1));
        }
    }

    return states_table;
}

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

DESystem op::SynchronizeStage2(op::StatesTableSTL const aTable, DESystem &aSys0,
                               DESystem &aSys1) {
    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    std::set<float> sync_events;
    std::set_union(aSys0.events_.begin(), aSys0.events_.end(),
                   aSys1.events_.begin(), aSys1.events_.end(),
                   std::inserter(sync_events, sync_events.begin()));

    DESystem::GraphHostData result(aTable.size(), aTable.size());

    for (auto q = aTable.begin(); q != aTable.end(); ++q) {
        for (auto event : sync_events) {
            bool const is_in_p =
                std::find(aSys0.events_.begin(), aSys0.events_.end(), event) !=
                aSys0.events_.end();
            bool const is_in_e =
                std::find(aSys1.events_.begin(), aSys1.events_.end(), event) !=
                aSys1.events_.end();
            if (!is_in_p && !is_in_e) {
                break;
            }

            int xid = -1;
            int yid;

            auto index0 = std::distance(aTable.begin(), q);

            if (is_in_p) {
                /*
                 * TODO: ublas::compressed_matrix<>::iterator1 does not have
                 * operators + and += implemented. Report bug and remove the
                 * following workaround.
                 */
                auto p_row = aSys0.graph_.begin1();
                for (auto i = 0; i < std::get<0>(*q); ++i) {
                    ++p_row;
                }
                // TODO: End of the workaround

                for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
                    if (std::fmod(*elem, event) == 0.0f) {
                        xid = elem.index2();
                        if (!is_in_e) {
                            yid = std::get<1>(*q);
                            StatesTupleSTL to_state_tuple =
                                std::make_tuple(xid, yid);
                            auto index1_iter = std::find_if(
                                aTable.begin(), aTable.end(),
                                [&to_state_tuple](
                                    std::tuple<cldes_size_t, cldes_size_t> i) {
                                    return (i == to_state_tuple);
                                });
                            size_t index1 =
                                std::distance(aTable.begin(), index1_iter);
                            if (result(index0, index1) != 0) {
                                result(index0, index1) *= event;
                            } else {
                                result(index0, index1) = event;
                            }
                            break;
                        }
                    }
                }
            }

            /*
             * TODO: ublas::compressed_matrix<>::iterator1 does not have
             * operators + and += implemented. Report bug and remove the
             * following workaround.
             */
            auto e_row = aSys1.graph_.begin1();
            for (auto i = 0; i < std::get<1>(*q); ++i) {
                ++e_row;
            }
            // TODO: End of the workaround

            if ((is_in_e && !is_in_p) || (is_in_e && (xid != -1))) {
                for (auto elem = e_row.begin(); elem != e_row.end(); ++elem) {
                    if (std::fmod(*elem, event) == 0.0f) {
                        if (!is_in_p) {
                            xid = std::get<0>(*q);
                        }
                        yid = elem.index2();
                        StatesTupleSTL to_state_tuple =
                            std::make_tuple(xid, yid);
                        auto index1_iter = std::find_if(
                            aTable.begin(), aTable.end(),
                            [&to_state_tuple](
                                std::tuple<cldes_size_t, cldes_size_t> i) {
                                return (i == to_state_tuple);
                            });
                        size_t index1 =
                            std::distance(aTable.begin(), index1_iter);
                        if (result(index0, index1) != 0) {
                            result(index0, index1) *= event;
                        } else {
                            result(index0, index1) = event;
                        }
                        break;
                    }
                }
            }
        }
    }

    // Copy device graph to host memory
    DESystem sync_sys(result, aTable.size(), initstate_sync, markedstates_sync);
    sync_sys.events_ = sync_events;

    return sync_sys;
}

op::StatesTupleSTL *op::TransitionVirtual(DESystem const &aP,
                                          DESystem const &aE,
                                          op::StatesTupleSTL const q,
                                          float const event) {
    bool const is_in_p = aP.events_.find(event) != aP.events_.end();
    bool const is_in_e = aE.events_.find(event) != aE.events_.end();

    if (!is_in_p && !is_in_e) {
        return nullptr;
    }

    /*
     * TODO: ublas::compressed_matrix<>::iterator1 does not have
     * operators + and += implemented. Report bug and remove the
     * following workaround.
     */
    auto p_row = aP.graph_.begin1();
    for (auto i = 0; i < std::get<0>(q); ++i) {
        ++p_row;
    }
    // TODO: End of the workaround

    StatesTupleSTL *ret = nullptr;

    if (is_in_p) {
        for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
            if (std::fmod(*elem, event) == 0.0f) {
                ret = new StatesTupleSTL;
                std::get<0>(*ret) = elem.index2();
                if (!is_in_e) {
                    std::get<1>(*ret) = std::get<1>(q);
                    return ret;
                }
            }
        }
    }

    /*
     * TODO: ublas::compressed_matrix<>::iterator1 does not have
     * operators + and += implemented. Report bug and remove the
     * following workaround.
     */
    auto e_row = aE.graph_.begin1();
    for (auto i = 0; i < std::get<1>(q); ++i) {
        ++e_row;
    }
    // TODO: End of the workaround

    if ((is_in_e && !is_in_p) || (is_in_e && (ret != nullptr))) {
        for (auto elem = e_row.begin(); elem != e_row.end(); ++elem) {
            if (std::fmod(*elem, event) == 0.0f) {
                if (!is_in_p) {
                    ret = new StatesTupleSTL;
                    std::get<0>(*ret) = std::get<0>(q);
                }
                std::get<1>(*ret) = elem.index2();
                return ret;
            }
        }
    }

    return ret;
}

bool op::TransitionReal(DESystem const &aP, cldes_size_t const &x,
                        float const &event) {
    bool const is_in_p = aP.events_.find(event) != aP.events_.end();

    if (!is_in_p) {
        return false;
    }

    /*
     * TODO: ublas::compressed_matrix<>::iterator1 does not have
     * operators + and += implemented. Report bug and remove the
     * following workaround.
     */
    auto p_row = aP.graph_.begin1();
    for (auto i = 0; i < x; ++i) {
        ++p_row;
    }
    // TODO: End of the workaround

    for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
        if (std::fmod(*elem, event) == 0.0f) {
            return true;
        }
    }

    return false;
}

template <class EventsType, class GraphType>
static op::StatesTableSTL __TransitionVirtualInv(EventsType const &aEventsP,
                                                 EventsType const &aEventsE,
                                                 GraphType const &aInvGraphP,
                                                 GraphType const &aInvGraphE,
                                                 op::StatesTupleSTL const q,
                                                 float const event) {
    bool const is_in_p = aEventsP.find(event) != aEventsP.end();
    bool const is_in_e = aEventsE.find(event) != aEventsE.end();

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
    for (auto i = 0; i < std::get<0>(q); ++i) {
        ++p_row;
    }
    auto e_row = aInvGraphE.begin1();
    for (auto i = 0; i < std::get<1>(q); ++i) {
        ++e_row;
    }
    // TODO: End of the workaround

    if (is_in_p) {
        for (auto elem = p_row.begin(); elem != p_row.end(); ++elem) {
            if (std::fmod(*elem, event) == 0.0f) {
                if (!is_in_e) {
                    ret.push_back(
                        std::make_tuple(elem.index2(), std::get<1>(q)));
                    return ret;
                } else {
                    for (auto elemg2 = e_row.begin(); elemg2 != e_row.end();
                         ++elemg2) {
                        if (std::fmod(*elemg2, event) == 0.0f) {
                            ret.push_back(std::make_tuple(elem.index2(),
                                                          elemg2.index2()));
                        }
                    }
                }
            }
        }
    } else {  // Is only on e: !is_in_p && is_in_e
        for (auto elemg2 = e_row.begin(); elemg2 != e_row.end(); ++e_row) {
            if (std::fmod(*elemg2, event) == 0.0f) {
                ret.push_back(std::make_tuple(std::get<0>(q), elemg2.index2()));
                return ret;
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
                              std::set<op::StatesTupleSTL> &C,
                              op::StatesTupleSTL const &q,
                              std::set<float> const &s_non_contr) {
    std::set<op::StatesTupleSTL> f;
    f.emplace(q);

    while (f.size() != 0) {
        op::StatesTupleSTL x = *(f.begin());
        f.erase(f.begin());
        C.erase(x);
        for (auto e : s_non_contr) {
            auto finv = __TransitionVirtualInv(aEventsP, aEventsE, aInvGraphP,
                                               aInvGraphE, q, e);
            if (finv.size() != 0) {
                f.insert(finv.begin(), finv.end());
            }
        }
    }

    return;
}

DESystem op::SupervisorSynth(DESystem &aP, DESystem &aE,
                             std::set<float> const &non_contr) {
    std::set<StatesTupleSTL> s_states;
    std::set<StatesTupleSTL> f;

    f.insert(std::make_tuple(aP.init_state_, aE.init_state_));

    DESystem::EventsSet s_events;
    std::set_union(aP.events_.begin(), aP.events_.end(), aE.events_.begin(),
                   aE.events_.end(), std::inserter(s_events, s_events.begin()));

    DESystem::EventsSet s_non_contr;
    std::set_intersection(s_events.begin(), s_events.end(), non_contr.begin(),
                          non_contr.end(),
                          std::inserter(s_non_contr, s_non_contr.begin()));

    auto const p_inversedgraph = boost::numeric::ublas::trans(aP.graph_);
    auto const e_inversedgraph = boost::numeric::ublas::trans(aE.graph_);

    while (f.size() != 0) {
        auto q = *(f.begin());
        s_states.emplace(q);
        f.erase(f.begin());
        for (auto event : s_events) {
            bool const is_non_contr =
                s_non_contr.find(event) != s_non_contr.end();
            auto fs_qevent = TransitionVirtual(aP, aE, q, event);

            if (is_non_contr && fs_qevent == nullptr &&
                TransitionReal(aP, std::get<0>(q), event)) {
                __RemoveBadStates(aP.events_, aE.events_, p_inversedgraph,
                                  e_inversedgraph, s_states, q, s_non_contr);
            } else if (fs_qevent) {
                auto is_in_f = f.find(*fs_qevent) != f.end();
                auto is_in_s_states =
                    s_states.find(*fs_qevent) != s_states.end();
                if (!is_in_s_states && !is_in_f) {
                    f.emplace(*fs_qevent);
                }
            }

            if (fs_qevent != nullptr) {
                delete fs_qevent;
            }
        }
    }

    StatesTableSTL s_virtual;
    std::copy(s_states.begin(), s_states.end(), std::back_inserter(s_virtual));
    auto supervisor = SynchronizeStage2(s_virtual, aP, aE);

    return supervisor.Trim();
}
