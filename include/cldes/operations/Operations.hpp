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

#include "cldes/constants.hpp"
#include <Eigen/Sparse>
#include <stack>
#include <tuple>

namespace cldes {

namespace op {

template<size_t NEvents, typename StorageIndex>
using GraphType =
  Eigen::SparseMatrix<std::bitset<NEvents>, Eigen::RowMajor, StorageIndex>;

template<typename StorageIndex>
using StatesArray = std::vector<StorageIndex>;

/*! \brief tuple representing a state of a virtual synch (stage 1)
 *
 * (state_id_g0, state_id_g1)
 */
template<typename StorageIndex>
using StatesTupleSTL = std::pair<StorageIndex, StorageIndex>;

/*! \brief Hash table of tuples representing a virtual synch (stage 1)
 */
template<typename StorageIndex>
using StatesTableSTL = spp::sparse_hash_set<StorageIndex>;

template<typename StorageIndex>
using SparseStatesMap = spp::sparse_hash_map<StorageIndex, StorageIndex>;

template<typename StorageIndex>
using StatesStack = std::stack<StorageIndex>;

using EventsTableSTL = spp::sparse_hash_set<unsigned>;
/*! \brief Returns the parallel composition between aSys0 and aSys1
 */
template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
Synchronize(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
            cldes::DESystem<NEvents, StorageIndex> const& aSys1);

template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SynchronizeStage1(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

template<size_t NEvents, typename StorageIndex>
void
SynchronizeStage2(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

template<size_t NEvents, typename StorageIndex>
StorageIndex
TransitionVirtual(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1,
                  StorageIndex const& q,
                  cldes::ScalarType const& event);

template<size_t NEvents, typename StorageIndex>
void
RemoveBadStates(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aE,
                GraphType<NEvents, StorageIndex> const& aInvGraphP,
                GraphType<NEvents, StorageIndex> const& aInvGraphE,
                StatesTableSTL<StorageIndex>& C,
                StorageIndex const& q,
                std::bitset<NEvents> const& bit_non_contr,
                StatesTableSTL<StorageIndex>& rmtable);

template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SupervisorSynth(cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aS,
                EventsTableSTL const& non_contr);

/*
cldes::DESystem<NEvents>CL SynchronizeStage2(StatesTable const *aTable,
                                    cldes::DESystem<NEvents>CL &aSys0,
                                    cldes::DESystem<NEvents>CL &aSys1);

StatesTable *SynchronizeStage1(cldes::DESystem<NEvents>CL const &aSys0,
                               cldes::DESystem<NEvents>CL const &aSys1);
typedef struct StatesTuple {
    cl_uint x0;
    cl_uint x1;
} StatesTuple;

typedef struct StatesTable {
    cl_uint tsize;
    StatesTuple *table;
} StatesTable;

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

cldes::DESystem<NEvents>CL Synchronize(cldes::DESystem<NEvents>CL &aSys0,
                              cldes::DESystem<NEvents>CL &aSys1);
*/


} // namespace op
} // namespace cldes

// Including implementation
#include "cldes/src/operations/OperationsCore.hpp"
#endif // DESYSTEM_HPP
