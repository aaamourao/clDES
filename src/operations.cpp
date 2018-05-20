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
#include <set>
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
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Sync stage 1
    auto states_tuple = new StatesTuple[table_size];
    for (unsigned int ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
        for (unsigned int ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
            unsigned int index = ix1 * aSys0.graph_->size1() + ix0;

            states_tuple[index].x0 = ix0;
            states_tuple[index].x1 = ix1;
        }
    }

    // Sync stage 2
    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    auto sys0_events = CalcEventsInt_(aSys0.events_);
    auto sys1_events = CalcEventsInt_(aSys1.events_);

    auto gcd_private = CalcGCD_(sys0_events, sys1_events);
    auto sys0_private = sys0_events / gcd_private;
    auto sys1_private = sys1_events / gcd_private;

    DESystem::GraphHostData result(table_size, table_size);

    for (auto i = 0; i < table_size; ++i) {
        auto state = states_tuple[i];

        /*
         * TODO: ublas::compressed_matrix<>::iterator1 does not have operators +
         * and += implemented. Report bug and remove the following workaround.
         */
        auto sys0_row = aSys0.graph_->begin1();
        for (auto i_sys0 = 0; i_sys0 != state.x0; ++i_sys0) {
            ++sys0_row;
        }
        // TODO: End of the workaround

        for (auto j = sys0_row.begin(); j != sys0_row.end(); ++j) {
            float sys0_elem = *j;
            if (sys0_private > 1.0f && sys0_elem > 1.0f) {
                float sys0_gcd_priv = CalcGCD_(sys0_private, sys0_elem);
                if (sys0_gcd_priv > 1.0f) {
                    unsigned int index0 =
                        state.x1 * aSys0.graph_->size1() + j.index2();
                    if (result(i, index0) > 1.0f) {
                        result(i, index0) *= sys0_gcd_priv;
                    } else {
                        result(i, index0) = sys0_gcd_priv;
                    }
                }
            }
            if (sys0_elem > 1.0f) {
                /*
                 * TODO: ublas::compressed_matrix<>::iterator1 does not have
                 * operators + and += implemented. Report bug and remove the
                 * following workaround.
                 */
                auto sys1_row = aSys1.graph_->begin1();
                for (auto i_sys1 = 0; i_sys1 != state.x1; ++i_sys1) {
                    ++sys1_row;
                }
                // TODO: End of workaround

                for (auto j_sys1 = sys1_row.begin(); j_sys1 != sys1_row.end();
                     ++j_sys1) {
                    float sys1_elem = *j_sys1;
                    if (sys1_private > 1.0f) {
                        float sys1_gcd_priv = CalcGCD_(sys1_private, sys1_elem);
                        if (sys1_gcd_priv > 1.0f) {
                            unsigned int index0 =
                                j_sys1.index2() * aSys0.graph_->size1() +
                                state.x0;
                            if (result(i, index0) > 1.0f) {
                                result(i, index0) *= sys1_gcd_priv;
                            } else {
                                result(i, index0) = sys1_gcd_priv;
                            }
                        }
                    }
                    float sync_gcd = CalcGCD_(sys0_elem, sys1_elem);
                    if (sync_gcd > 1.0f) {
                        unsigned int index0 =
                            j_sys1.index2() * aSys0.graph_->size1() +
                            j.index2();
                        if (result(i, index0) > 1.0f) {
                            result(i, index0) *= sync_gcd;
                        } else {
                            result(i, index0) = sync_gcd;
                        }
                    }
                }
            }
        }
    }

    // Copy device graph to host memory
    DESystem sync_sys(result, table_size, initstate_sync, markedstates_sync);

    return sync_sys;
}

op::StatesTable *op::SynchronizeStage1(DESystem const &aSys0,
                                       DESystem const &aSys1) {
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Get the result and saves it on host memory
    auto states_tuple = new StatesTuple[table_size];
    auto sys0_size = aSys0.graph_->size1();

    for (auto ix0 = 0; ix0 < aSys0.states_number_; ++ix0) {
        for (auto ix1 = 0; ix1 < aSys1.states_number_; ++ix1) {
            auto index = ix1 * sys0_size + ix0;

            states_tuple[index].x0 = ix0;
            states_tuple[index].x1 = ix1;
        }
    }

    auto states_table = new StatesTable;
    states_table->tsize = table_size;
    states_table->table = states_tuple;

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

DESystem op::SynchronizeStage2(op::StatesTable const *aTable, DESystem &aSys0,
                               DESystem &aSys1) {
    auto initstate_sync = op::TablePos_(aSys0.init_state_, aSys1.init_state_,
                                        aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    auto sys0_events = CalcEventsInt_(aSys0.events_);
    auto sys1_events = CalcEventsInt_(aSys1.events_);

    auto gcd_private = CalcGCD_(sys0_events, sys1_events);
    auto sys0_private = sys0_events / gcd_private;
    auto sys1_private = sys1_events / gcd_private;

    DESystem::GraphHostData result(aTable->tsize, aTable->tsize);

    for (auto i = 0; i < aTable->tsize; ++i) {
        auto state = aTable->table[i];

        /*
         * TODO: ublas::compressed_matrix<>::iterator1 does not have operators +
         * and += implemented. Report bug and remove the following workaround.
         */
        auto sys0_row = aSys0.graph_->begin1();
        for (auto i_sys0 = 0; i_sys0 != state.x0; ++i_sys0) {
            ++sys0_row;
        }
        // TODO: End of the workaround

        for (auto j = sys0_row.begin(); j != sys0_row.end(); ++j) {
            float sys0_elem = *j;
            if (sys0_private > 1.0f && sys0_elem > 1.0f) {
                float sys0_gcd_priv = CalcGCD_(sys0_private, sys0_elem);
                if (sys0_gcd_priv > 1.0f) {
                    unsigned int index0 =
                        state.x1 * aSys0.graph_->size1() + j.index2();
                    if (result(i, index0) > 1.0f) {
                        result(i, index0) *= sys0_gcd_priv;
                    } else {
                        result(i, index0) = sys0_gcd_priv;
                    }
                }
            }
            if (sys0_elem > 1.0f) {
                /*
                 * TODO: ublas::compressed_matrix<>::iterator1 does not have
                 * operators + and += implemented. Report bug and remove the
                 * following workaround.
                 */
                auto sys1_row = aSys1.graph_->begin1();
                for (auto i_sys1 = 0; i_sys1 != state.x1; ++i_sys1) {
                    ++sys1_row;
                }
                // TODO: End of the workaround

                for (auto j_sys1 = sys1_row.begin(); j_sys1 != sys1_row.end();
                     ++j_sys1) {
                    float sys1_elem = *j_sys1;
                    if (sys1_private > 1.0f) {
                        float sys1_gcd_priv = CalcGCD_(sys1_private, sys1_elem);
                        if (sys1_gcd_priv > 1.0f) {
                            unsigned int index0 =
                                j_sys1.index2() * aSys0.graph_->size1() +
                                state.x0;
                            if (result(i, index0) > 1.0f) {
                                result(i, index0) *= sys1_gcd_priv;
                            } else {
                                result(i, index0) = sys1_gcd_priv;
                            }
                        }
                    }
                    float sync_gcd = CalcGCD_(sys0_elem, sys1_elem);
                    if (sync_gcd > 1.0f) {
                        unsigned int index0 =
                            j_sys1.index2() * aSys0.graph_->size1() +
                            j.index2();
                        if (result(i, index0) > 1.0f) {
                            result(i, index0) *= sync_gcd;
                        } else {
                            result(i, index0) = sync_gcd;
                        }
                    }
                }
            }
        }
    }

    // Copy device graph to host memory
    DESystem sync_sys(result, aTable->tsize, initstate_sync, markedstates_sync);

    return sync_sys;
}
