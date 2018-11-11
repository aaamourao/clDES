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

 LacSED - Laboratorio de Analise e Controle de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: cldes/src/des/DESystemCLCore.hpp
 Description: DESystemCL class definition.
 =========================================================================
*/

#include <algorithm>
#include <functional>
#include <vector>

namespace cldes {

// OclBackend disable by now
template<uint8_t NEvents, typename StorageIndex>
std::unique_ptr<backend::OclBackend>
  DESystemCL<NEvents, StorageIndex>::backend_ptr_ = nullptr;

template<uint8_t NEvents, typename StorageIndex>
DESystemCL<NEvents, StorageIndex>::DESystemCL(
  DESystem<NEvents, StorageIndex>& aSys)
  : DESystemBase{ aSys.Size(), aSys.GetInitialState() }
{
    GraphHostData ident{ this->states_number_, this->states_number_ };
    ident.setIdentity();
    ident += graph_;
    graph_ = ident.template cast<float>();
    this->marked_states_ = aSys.GetMarkedStates();
}

template<uint8_t NEvents, typename StorageIndex>
std::shared_ptr<DESystemBase<NEvents, StorageIndex>>
DESystemCL<NEvents, StorageIndex>::Clone() const
{
    std::shared_ptr<DESystemBase> this_ptr =
      std::make_shared<DESystemCL>(*this);
    return this_ptr;
}

template<uint8_t NEvents, typename StorageIndex>
bool
DESystemCL<NEvents, StorageIndex>::IsVirtual() const
{
    return false;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystemCL<NEvents, StorageIndex>::GraphDeviceData
DESystemCL<NEvents, StorageIndex>::GetGraph() const
{
    return graph_;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystemCL<NEvents, StorageIndex>::StatesSet
DESystemCL<NEvents, StorageIndex>::accessiblePart()
{
    auto accessible_states = bfs_();
    return *accessible_states;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystemCL<NEvents, StorageIndex>::StatesSet
DESystemCL<NEvents, StorageIndex>::coaccessiblePart()
{
    StorageIndexSigned const n_marked = this->marked_states_.size();
    StatesVector host_x{ static_cast<StorageIndexSigned>(this->states_number_),
                         n_marked };
    host_x.reserve(this->marked_states_.size());

    {
        auto pos = 0ul;
        for (auto state : this->marked_states_) {
            host_x.coeffRef(state, pos) = true;
            ++pos;
        }
    }
    // Copy search vector to device memory
    StatesDeviceVector x;
    viennacl::copy(host_x, x);

    // Copy to device memory
    viennacl::copy(graph_.transpose().eval(), device_graph_);

    // Executes BFS
    StatesDeviceVector y;
    auto n_accessed_states = 0ul;
    for (auto i = 0ul; i < this->states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = viennacl::linalg::prod(device_graph_, x);

        x = std::move(y);

        if (n_accessed_states == x.nnz()) {
            break;
        } else {
            n_accessed_states = x.nnz();
        }
    }

    viennacl::copy(x, host_x);

    auto coaccessible_states = StatesSet{};

    // Add results to a std::set vector
    for (StorageIndexSigned s = 0; s < host_x.outerSize(); ++s) {
        for (ColIteratorConst e(host_x, s); e; ++e) {
            coaccessible_states.emplace(e.row());
        }
    }

    return coaccessible_states;
}

template<uint8_t NEvents, typename StorageIndex>
template<class StatesType>
std::shared_ptr<typename DESystemCL<NEvents, StorageIndex>::StatesSet>
DESystemCL<NEvents, StorageIndex>::bfs_(
  StatesType const& aInitialNode,
  std::function<void(StorageIndex const&, StorageIndex const&)> const&
    aBfsVisit)
{
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    StatesVector host_x{ static_cast<Eigen::Index>(this->states_number_), 1 };

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x.coeffRef(aInitialNode, 0) = 1;

    return bfsCalc_(host_x, aBfsVisit, nullptr);
}

template<uint8_t NEvents, typename StorageIndex>
std::shared_ptr<typename DESystemCL<NEvents, StorageIndex>::StatesSet>
DESystemCL<NEvents, StorageIndex>::bfsCalc_(
  StatesVector& aHostX,
  std::function<void(StorageIndex const&, StorageIndex const&)> const&
    aBfsVisit,
  std::vector<StorageIndex> const* const aStatesMap)
{
    cl_uint n_initial_nodes = aHostX.cols();

    // Copy search vector to device memory
    StatesDeviceVector x;
    viennacl::copy(aHostX, x);

    // Copy to device memory
    viennacl::copy(graph_, device_graph_);

    // Executes BFS
    StatesDeviceVector y;
    auto n_accessed_states = 0ul;
    for (auto i = 0u; i < this->states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = viennacl::linalg::prod(device_graph_, x);

        x = std::move(y);

        if (n_accessed_states == x.nnz()) {
            break;
        } else {
            n_accessed_states = x.nnz();
        }
    }

    viennacl::copy(x, aHostX);

    // Unfortunatelly, only C++17 allows shared_ptr to arrays
    std::shared_ptr<StatesSet> accessed_states{
        new StatesSet[n_initial_nodes], std::default_delete<StatesSet[]>()
    };

    // Add results to a std::set vector
    for (StorageIndexSigned s = 0; s < aHostX.outerSize(); ++s) {
        for (ColIteratorConst e(aHostX, s); e; ++e) {
            if (!aStatesMap && !aBfsVisit) {
                accessed_states.get()[e.col()].emplace(e.row());
            }
        }
    }

    return accessed_states;
}

template<uint8_t NEvents, typename StorageIndex>
std::shared_ptr<typename DESystemCL<NEvents, StorageIndex>::StatesSet>
DESystemCL<NEvents, StorageIndex>::bfs_()
{
    return bfs_(this->init_state_, nullptr);
}

template<uint8_t NEvents, typename StorageIndex>
bool
DESystemCL<NEvents, StorageIndex>::Containstrans(StorageIndex const&,
                                                 ScalarType const&) const
{
    return false;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystemCL<NEvents, StorageIndex>::StorageIndexSigned
DESystemCL<NEvents, StorageIndex>::trans(StorageIndex const&,
                                         ScalarType const&) const
{
    return 0;
}

template<uint8_t NEvents, typename StorageIndex>
bool
DESystemCL<NEvents, StorageIndex>::Containsinvtrans(StorageIndex const&,
                                                    ScalarType const&) const
{
    return false;
}

template<uint8_t NEvents, typename StorageIndex>
StatesArray<StorageIndex>
DESystemCL<NEvents, StorageIndex>::invtrans(StorageIndex const&,
                                            ScalarType const&) const
{
    StatesArray<StorageIndex> foo{};
    return foo;
}

template<uint8_t NEvents, typename StorageIndex>
EventsSet<NEvents>
DESystemCL<NEvents, StorageIndex>::GetStateEvents(StorageIndex const&) const
{
    EventsSet<NEvents> foo;
    return foo;
}

template<uint8_t NEvents, typename StorageIndex>
EventsSet<NEvents>
DESystemCL<NEvents, StorageIndex>::GetInvStateEvents(StorageIndex const&) const
{
    EventsSet<NEvents> foo;
    return foo;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemCL<NEvents, StorageIndex>::AllocateInvertedGraph() const
{
    return;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemCL<NEvents, StorageIndex>::ClearInvertedGraph() const
{
    return;
}
}
