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
#include <CL/cl.hpp>
#include <algorithm>
#include <set>
#include "backend/oclbackend.hpp"
#include "des/desystem.hpp"

using namespace cldes;

cldes::ScalarType op::GetTransitions(DESystem const &aSys) {
    auto graph = aSys.GetGraph();
    float transitions = 1.0f;

    for (auto lin = graph.begin1(); lin != graph.end1(); ++lin) {
        for (auto col = lin.begin(); col != lin.end(); ++col) {
            if (int(transitions) % int((*col))) {
                transitions = transitions * (*col);
            }
        }
    }

    return transitions;
}

ScalarType op::GetPrivateTransitions(DESystem const &aSysTarget,
                                     DESystem const &aSys) {
    auto transitionsTarget = op::GetTransitions(aSysTarget);
    auto graph = aSys.GetGraph();

    auto privatetransitions = 1.0f;

    for (auto lin = graph.begin1(); lin != graph.end1(); ++lin) {
        for (auto col = lin.begin(); col != lin.end(); ++col) {
            if (int(transitionsTarget) % int(*col) &&
                int(privatetransitions) % int(*col)) {
                privatetransitions = privatetransitions * (*col);
            }
        }
    }

    auto transitions = op::GetTransitions(aSys);
    auto graphTarget = aSysTarget.GetGraph();

    for (auto lin = graphTarget.begin1(); lin != graphTarget.end1(); ++lin) {
        for (auto col = lin.begin(); col != lin.end(); ++col) {
            if (int(transitions) % int(*col) &&
                int(privatetransitions) % int(*col)) {
                privatetransitions = privatetransitions * (*col);
            }
        }
    }

    return privatetransitions;
}

DESystem op::Synchronize(DESystem const &aSys0, DESystem const &aSys1) {
    auto oclbackend = backend::OclBackend::Instance();

    auto new_size = aSys0.states_number_ * aSys1.states_number_;

    // Get new marked states
    std::set<cldes_size_t> marked_states;
    std::set_union(aSys0.marked_states_.begin(), aSys0.marked_states_.end(),
                   aSys1.marked_states_.begin(), aSys1.marked_states_.end(),
                   std::inserter(marked_states, marked_states.begin()));

    // Get new initial states
    auto init_state = aSys0.init_state_;

    // Create new system with empty transitions
    DESystem compdes{new_size, 0, marked_states};

    // TODO: Calculate the transitions

    return aSys1;
}

op::StatesTuple *op::SynchronizeStage1(DESystem const &aSys0,
                                       DESystem const &aSys1) {
    auto oclbackend = backend::OclBackend::Instance();

    auto new_size = aSys0.states_number_ * aSys1.states_number_;

    auto states_tuple_dev = oclbackend->CreateBuffer(
        CL_MEM_WRITE_ONLY, new_size * sizeof(StatesTuple), nullptr);

    auto syncstage1kernel = oclbackend->GetKernel("Synchronize_Stage1");

    syncstage1kernel.local_work_size(0, 1);
    syncstage1kernel.local_work_size(1, 1);
    syncstage1kernel.global_work_size(0, aSys0.states_number_);
    syncstage1kernel.global_work_size(1, aSys1.states_number_);

    oclbackend->Enqueue(syncstage1kernel(
        states_tuple_dev, static_cast<cl_uint>(aSys0.states_number_)));

    auto states_tuple_host = new StatesTuple[new_size];
    clEnqueueReadBuffer(oclbackend->CommandQueue(), states_tuple_dev, CL_TRUE,
                        0, new_size * sizeof(StatesTuple),
                        (void *)states_tuple_host, 0, NULL, NULL);

    return states_tuple_host;
}
