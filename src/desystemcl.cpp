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

 File: desystem.cpp
 Description: DESystemCL class implementation. DESystemCL is a graph, which
 is modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <functional>
#include <vector>
#include "backend/oclbackend.hpp"
#include "des/desystemcl.hpp"
#include "operations/operations.hpp"
#include "viennacl/linalg/prod.hpp"

using namespace cldes;

DESystemCL::DESystemCL(GraphHostData const &aGraph,
                       cldes_size_t const &aStatesNumber,
                       cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                       bool const &aDevCacheEnabled)
    : graph_{new GraphHostData{aGraph}}, init_state_{aInitState} {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;
    device_graph_ = nullptr;

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        this->CacheGraph_();
    }
}

DESystemCL::DESystemCL(cldes_size_t const &aStatesNumber,
                       cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                       bool const &aDevCacheEnabled)
    : DESystemCL::DESystemCL{GraphHostData{aStatesNumber, aStatesNumber},
                             aStatesNumber, aInitState, aMarkedStates,
                             aDevCacheEnabled} {}

DESystemCL::DESystemCL(DESystemCL const &aSys)
    : graph_{new GraphHostData{*(aSys.graph_)}}, init_state_{aSys.init_state_} {
    states_number_ = aSys.states_number_;
    marked_states_ = aSys.marked_states_;
    dev_cache_enabled_ = aSys.dev_cache_enabled_;
    is_cache_outdated_ = aSys.is_cache_outdated_;

    if (device_graph_) {
        device_graph_ = new GraphDeviceData{*(aSys.device_graph_)};
    } else {
        device_graph_ = nullptr;
    }
}

DESystemCL::~DESystemCL() {
    // Delete uBlas data
    if (graph_) {
        delete graph_;
    }
    if (dev_cache_enabled_) {
        delete device_graph_;
    }
}

DESystemCL::GraphHostData DESystemCL::GetGraph() const { return *graph_; }

void DESystemCL::CacheGraph_() {
    // Allocate space for device_graph_
    if (device_graph_ == nullptr) {
        device_graph_ = new viennacl::compressed_matrix<ScalarType>{
            states_number_, states_number_};
    }
    viennacl::copy(trans(*graph_), *device_graph_);

    is_cache_outdated_ = false;
}

void DESystemCL::UpdateGraphCache_() {
    viennacl::copy(trans(*graph_), *device_graph_);
    is_cache_outdated_ = false;
}

template <typename StatesType>
DESystemCL::StatesSet *DESystemCL::Bfs_(
    StatesType const &aInitialNodes,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit) {
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    // There is no need of search if a marked state is coaccessible
    StatesVector host_x{states_number_, aInitialNodes.size()};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    std::vector<cldes_size_t> states_map;
    for (auto state : aInitialNodes) {
        host_x(state, states_map.size()) = 1;
        // Maping each search from each initial node to their correspondent
        // vector on the matrix
        states_map.push_back(state);
    }

    return BfsCalc_(host_x, aBfsVisit, &states_map);
}

template <>
DESystemCL::StatesSet *DESystemCL::Bfs_<cldes_size_t>(
    cldes_size_t const &aInitialNode,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit) {
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    StatesVector host_x{states_number_, 1};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x(aInitialNode, 0) = 1;

    return BfsCalc_(host_x, aBfsVisit, nullptr);
}

DESystemCL::StatesSet *DESystemCL::Bfs_() {
    return Bfs_(init_state_, nullptr);
};

DESystemCL::StatesSet *DESystemCL::BfsCalc_(
    StatesVector &aHostX,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit,
    std::vector<cldes_size_t> const *const aStatesMap) {
    cl_uint n_initial_nodes = aHostX.size2();

    // Copy search vector to device memory
    StatesDeviceVector x;
    viennacl::copy(aHostX, x);

    ublas::compressed_matrix<ScalarType> graph{*graph_};

    // Sum it to identity and set all nnz to 1.0f
    for (auto i = graph.begin1(); i != graph.end1(); ++i) {
        graph(i.index1(), i.index2()) = 1.0f;
    }

    // Copy to device memory
    GraphDeviceData dev_graph;
    float pattern = 1.0f;
    viennacl::copy(trans(graph), dev_graph);
    clEnqueueFillBuffer(viennacl::ocl::get_queue().handle().get(),
                        dev_graph.handle().opencl_handle(), &pattern,
                        sizeof(float), 0, dev_graph.nnz() * sizeof(float), 0,
                        NULL, NULL);

    // Executes BFS
    StatesDeviceVector y;
    auto n_accessed_states = 0;
    for (auto i = 0; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = viennacl::linalg::prod(dev_graph, x);

        if (n_accessed_states == y.nnz()) {
            break;
        } else {
            n_accessed_states = y.nnz();
        }

        x = y;
    }
    viennacl::copy(y, aHostX);
    // Add results to a std::set vector
    auto accessed_states = new StatesSet[n_initial_nodes];
    for (auto node = aHostX.begin1(); node != aHostX.end1(); ++node) {
        for (auto elem = node.begin(); elem != node.end(); ++elem) {
            accessed_states[elem.index2()].emplace(node.index1());
        }
    }
    return accessed_states;
}

DESystemCL::StatesSet *DESystemCL::BfsSpMV_(
    cldes_size_t const &aInitialNode,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit) {
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    StatesDenseVector host_x{states_number_};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x.clear();
    host_x[aInitialNode] = static_cast<ScalarType>(1);

    return BfsCalcSpMV_(host_x, aBfsVisit);
}

DESystemCL::StatesSet *DESystemCL::BfsCalcSpMV_(
    StatesDenseVector &aHostX,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit) {
    // Cache graph temporally
    if (!dev_cache_enabled_) {
        this->CacheGraph_();
    } else if (is_cache_outdated_) {
        this->UpdateGraphCache_();
    }

    // Copy search vector to device memory
    StatesDeviceDenseVector x{states_number_};
    viennacl::copy(aHostX, x);

    auto accessed_states = new StatesSet;

    // Executes BFS
    StatesDeviceDenseVector y;
    for (auto i = 0; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = viennacl::linalg::prod(*device_graph_, x);
        x = y;

        bool accessed_new_state = false;

        // Unfortunatelly, until now, ViennaCL does not allow iterating on
        // compressed matrices. Until it is implemented, it is necessary
        // to copy the vector to the host memory.
        viennacl::copy(x, aHostX);

        for (auto vert = aHostX.begin(); vert != aHostX.end(); ++vert) {
            // TODO: Use FUNCTIONAL FILTERING here
            if (*vert) {
                auto accessed_state = vert.index();
                if (accessed_states->emplace(accessed_state).second) {
                    accessed_new_state = true;
                    if (aBfsVisit) {
                        aBfsVisit(0, accessed_state);
                    }
                }
            }
        }

        // TODO: If some search did not accessed a new state, it does not need
        // to be in the BFS Matrix

        // If all accessed states were previously "colored", stop searching.
        if (!accessed_new_state) {
            break;
        }
    }

    // Remove graph_ from device memory, if it is set so
    if (!dev_cache_enabled_) {
        delete device_graph_;
    }

    return accessed_states;
}

DESystemCL::StatesSet DESystemCL::AccessiblePart() {
    // Executes a BFS on graph_
    auto paccessible_states = Bfs_();

    auto accessible_states = *paccessible_states;
    delete[] paccessible_states;

    return accessible_states;
}

DESystemCL::StatesSet DESystemCL::CoaccessiblePart() {
    StatesSet searching_nodes;
    // Initialize initial_nodes with all states, but marked states
    {
        StatesSet all_nodes(
            boost::counting_iterator<cldes_size_t>(0),
            boost::counting_iterator<cldes_size_t>(states_number_));
        std::set_difference(
            all_nodes.begin(), all_nodes.end(), marked_states_.begin(),
            marked_states_.end(),
            std::inserter(searching_nodes, searching_nodes.begin()));
    }

    StatesSet coaccessible_states = marked_states_;
    auto paccessed_states = Bfs_(
        searching_nodes,
        [this, &coaccessible_states](cldes_size_t const &aInitialState,
                                     cldes_size_t const &aAccessedState) {
            bool const is_marked = this->marked_states_.find(aAccessedState) !=
                                   this->marked_states_.end();
            if (is_marked) {
                coaccessible_states.emplace(aInitialState);
            }
            // TODO: If all states are marked, the search is done
        });
    delete[] paccessed_states;

    return coaccessible_states;
}

DESystemCL::StatesSet DESystemCL::TrimStates() {
    auto accpart = this->AccessiblePart();
    auto coaccpart = this->CoaccessiblePart();

    StatesSet trimstates;
    std::set_intersection(accpart.begin(), accpart.end(), coaccpart.begin(),
                          coaccpart.end(),
                          std::inserter(trimstates, trimstates.begin()));

    return trimstates;
}

DESystemCL DESystemCL::Trim(bool const &aDevCacheEnabled) {
    auto trimstates = this->TrimStates();

    if (trimstates.size() == graph_->size1()) {
        return *this;
    }

    // First remove rows of non-trim states
    GraphHostData sliced_graph{trimstates.size(), states_number_};
    auto mapped_state = 0;
    for (auto state : trimstates) {
        ublas::matrix_row<GraphHostData> graphrow(*graph_, state);
        ublas::matrix_row<GraphHostData> slicedrow(sliced_graph, mapped_state);
        slicedrow.swap(graphrow);
        ++mapped_state;
    }

    // Now remove the non-trim columns
    GraphHostData trim_graph{trimstates.size(), trimstates.size()};
    mapped_state = 0;
    for (auto state : trimstates) {
        ublas::matrix_column<GraphHostData> slicedcol(sliced_graph, state);
        ublas::matrix_column<GraphHostData> trimcol(trim_graph, mapped_state);
        trimcol.swap(slicedcol);
        ++mapped_state;
    }

    StatesSet trim_marked_states;
    std::set_intersection(
        marked_states_.begin(), marked_states_.end(), trimstates.begin(),
        trimstates.end(),
        std::inserter(trim_marked_states, trim_marked_states.begin()));

    DESystemCL trim_system{trim_graph, trimstates.size(), init_state_,
                           trim_marked_states, aDevCacheEnabled};
    return trim_system;
}

void DESystemCL::InsertEvents(DESystemCL::EventsSet const &aEvents) {
    events_ = EventsSet(aEvents);
}
