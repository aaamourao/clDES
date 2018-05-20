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

 File: operations.hpp
 Description: Declaration of operation functions.
 =========================================================================
*/

#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "constants.hpp"

namespace cldes {

class DESystemCL;
class DESystem;

namespace op {

typedef struct StatesTuple {
    cl_uint x0;
    cl_uint x1;
} StatesTuple;

typedef struct StatesTable {
    cl_uint tsize;
    StatesTuple *table;
} StatesTable;

cldes::cldes_size_t TablePos_(cldes::cldes_size_t const &aG0Pos,
                              cldes::cldes_size_t const &aG1Pos,
                              cldes::cldes_size_t const &aG0NStates);

template <class KernelType>
void SetWorkGroups_(KernelType *k, cldes::cldes_size_t const aGws0,
                    cldes::cldes_size_t const aGws1,
                    cldes::cldes_size_t const aLws0,
                    cldes::cldes_size_t const aLws1) {
    k->local_work_size(0, aLws0);
    k->local_work_size(1, aLws1);
    k->global_work_size(0, aGws0);
    k->global_work_size(1, aGws1);
}

template <class EventsSetType>
float CalcEventsInt_(EventsSetType const &aEvents) {
    float events_integer = 1.0f;

    for (auto event : aEvents) {
        events_integer = events_integer * event;
    }

    return events_integer;
}

float CalcGCD_(float aG0, float aG1);

cldes::DESystemCL Synchronize(cldes::DESystemCL &aSys0,
                              cldes::DESystemCL &aSys1);

cldes::DESystem Synchronize(cldes::DESystem &aSys0, cldes::DESystem &aSys1);

StatesTable *SynchronizeStage1(cldes::DESystemCL const &aSys0,
                               cldes::DESystemCL const &aSys1);

StatesTable *SynchronizeStage1(cldes::DESystem const &aSys0,
                               cldes::DESystem const &aSys1);

cldes::DESystemCL SynchronizeStage2(StatesTable const *aTable,
                                    cldes::DESystemCL &aSys0,
                                    cldes::DESystemCL &aSys1);

cldes::DESystem SynchronizeStage2(StatesTable const *aTable,
                                  cldes::DESystem &aSys0,
                                  cldes::DESystem &aSys1);
}  // namespace op
}  // namespace cldes
#endif  // DESYSTEM_HPP
