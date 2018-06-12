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

#include "backend/OCLBackend.hpp"
#include "cldes/DESystemCL.hpp"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>

DESystemCL<NEvents>
op::Synchronize(DESystemCL<NEvents>& aSys0, DESystemCL<NEvents>& aSys1)
{
    /*
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
      CL_MEM_READ_WRITE, table_size * sizeof(StatesTuple), nullptr);

    auto syncstage1kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage1");

    // Set Work groups size
    SetWorkGroups_(
      &syncstage1kernel, aSys0.states_number_, aSys1.states_number_, 1, 1);

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage1kernel(
      states_tuple_dev, static_cast<cl_uint>(aSys0.states_number_)));

    // Get the result and saves it on host memory
    auto states_tuple_host = new StatesTuple[table_size];
    clEnqueueReadBuffer(aSys0.backend_ptr_->CommandQueue(),
                        states_tuple_dev,
                        CL_TRUE,
                        0,
                        table_size * sizeof(StatesTuple),
                        (void*)states_tuple_host,
                        0,
                        NULL,
                        NULL);

    if (aSys0.is_cache_outdated_) {
        aSys0.UpdateGraphCache_();
    }

    if (aSys1.is_cache_outdated_) {
        aSys1.UpdateGraphCache_();
    }

    auto initstate_sync =
      op::TablePos_(aSys0.init_state_, aSys1.init_state_, aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(),
                   aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(),
                   aSys1.marked_states_.end(),
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
    aSys0.backend_ptr_->Enqueue(
      syncstage2kernel(states_tuple_dev,
                       aSys0.device_graph_->handle1().opencl_handle(),
                       aSys0.device_graph_->handle2().opencl_handle(),
                       aSys0.device_graph_->handle().opencl_handle(),
                       static_cast<cl_uint>(aSys0.device_graph_->rows()),
                       asys0_private,
                       aSys1.device_graph_->handle1().opencl_handle(),
                       aSys1.device_graph_->handle2().opencl_handle(),
                       aSys1.device_graph_->handle().opencl_handle(),
                       asys1_private,
                       result_dev.handle().opencl_handle(),
                       static_cast<cl_uint>(result_dev.internal_rows())));

    // Copy device graph to host memory
    DESystem<NEvents> CL sync_sys(
      table_size, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));
    */
    DESystemCL sync_sys{};

    return sync_sys;
}

op::StatesTable*
op::SynchronizeStage1(DESystemCL<NEvents> const& aSys0,
                      DESystemCL<NEvents> const& aSys1)
{
    /*
    auto table_size = aSys0.states_number_ * aSys1.states_number_;

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
      CL_MEM_WRITE_ONLY, table_size * sizeof(StatesTuple), nullptr);

    auto syncstage1kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage1");

    // Set Work groups size
    SetWorkGroups_(
      &syncstage1kernel, aSys0.states_number_, aSys1.states_number_, 1, 1);

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(syncstage1kernel(
      states_tuple_dev, static_cast<cl_uint>(aSys0.states_number_)));

    // Get the result and saves it on host memory
    auto states_tuple_host = new StatesTuple[table_size];
    clEnqueueReadBuffer(aSys0.backend_ptr_->CommandQueue(),
                        states_tuple_dev,
                        CL_TRUE,
                        0,
                        table_size * sizeof(StatesTuple),
                        (void*)states_tuple_host,
                        0,
                        NULL,
                        NULL);

    auto states_table = new StatesTable;
    states_table->tsize = table_size;
    states_table->table = states_tuple_host;

    return states_table;
    */
    return nullptr;
}

DESystemCL<NEvents>
op::SynchronizeStage2(op::StatesTable const* aTable,
                      DESystemCL<NEvents>& aSys0,
                      DESystemCL<NEvents>& aSys1)
{
    /*
    if (aSys0.is_cache_outdated_) {
        aSys0.UpdateGraphCache_();
    }

    if (aSys1.is_cache_outdated_) {
        aSys1.UpdateGraphCache_();
    }

    auto initstate_sync =
      op::TablePos_(aSys0.init_state_, aSys1.init_state_, aSys0.states_number_);

    std::set<cldes_size_t> markedstates_sync;
    std::set_union(aSys0.marked_states_.begin(),
                   aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(),
                   aSys1.marked_states_.end(),
                   std::inserter(markedstates_sync, markedstates_sync.begin()));

    // Allocate memory on the device
    auto states_tuple_dev = aSys0.backend_ptr_->GetContext().create_memory(
      CL_MEM_READ_ONLY, aTable->tsize * sizeof(StatesTuple), aTable->table);

    auto syncstage2kernel = aSys0.backend_ptr_->GetKernel("Synchronize_Stage2");

    // Set Work groups size
    SetWorkGroups_(
      &syncstage2kernel, aTable->tsize, aSys0.events_.size(), 1, 1);

    auto asys0_events = CalcEventsInt_(aSys0.events_);
    auto asys1_events = CalcEventsInt_(aSys1.events_);

    auto gcd_private = CalcGCD_(asys0_events, asys1_events);
    auto asys0_private = asys0_events / gcd_private;
    auto asys1_private = asys1_events / gcd_private;

    viennacl::matrix<float> result_dev(aTable->tsize, aTable->tsize);
    result_dev.clear();

    // Execute kernel on the device
    aSys0.backend_ptr_->Enqueue(
      syncstage2kernel(states_tuple_dev,
                       aSys0.device_graph_->handle1().opencl_handle(),
                       aSys0.device_graph_->handle2().opencl_handle(),
                       aSys0.device_graph_->handle().opencl_handle(),
                       static_cast<cl_uint>(aSys0.device_graph_->rows()),
                       asys0_private,
                       aSys1.device_graph_->handle1().opencl_handle(),
                       aSys1.device_graph_->handle2().opencl_handle(),
                       aSys1.device_graph_->handle().opencl_handle(),
                       asys1_private,
                       result_dev.handle().opencl_handle(),
                       static_cast<cl_uint>(result_dev.internal_rows())));

    // Copy device graph to host memory
    DESystem<NEvents> CL sync_sys(
      aTable->tsize, initstate_sync, markedstates_sync);
    viennacl::copy(result_dev, *(sync_sys.graph_));
    viennacl::copy(trans(*(sync_sys.graph_)), *(sync_sys.device_graph_));
    */
    DESystemCL sync_sys{};

    return sync_sys;
}
