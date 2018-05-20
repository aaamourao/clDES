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
 Description: DESystem class implementation. DESystem is a graph, which
 is modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#include "des/desystem.hpp"
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <functional>
#include <vector>
#include "backend/oclbackend.hpp"
#include "des/transition_proxy.hpp"
#include "operations/operations.hpp"

using namespace cldes;

DESystem::DESystem(GraphHostData const &aGraph,
                   cldes_size_t const &aStatesNumber,
                   cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : graph_{new GraphHostData{aGraph}}, init_state_{aInitState} {
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        this->CacheGraph_();
    }
}

DESystem::DESystem(cldes_size_t const &aStatesNumber,
                   cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                   bool const &aDevCacheEnabled)
    : DESystem::DESystem{GraphHostData{aStatesNumber, aStatesNumber},
                         aStatesNumber, aInitState, aMarkedStates,
                         aDevCacheEnabled} {}

DESystem::DESystem(DESystem const &aSys)
    : graph_{new GraphHostData{*(aSys.graph_)}}, init_state_{aSys.init_state_} {
    states_number_ = aSys.states_number_;
    marked_states_ = aSys.marked_states_;
    dev_cache_enabled_ = aSys.dev_cache_enabled_;
    is_cache_outdated_ = aSys.is_cache_outdated_;
}

DESystem::~DESystem() {
    // Delete uBlas data
    if (graph_) {
        delete graph_;
    }
}

DESystem::GraphHostData DESystem::GetGraph() const { return *graph_; }

void DESystem::CacheGraph_() { is_cache_outdated_ = false; }

DESystem::GraphHostData::const_reference DESystem::operator()(
    cldes_size_t const &aLin, cldes_size_t const &aCol) const {
    return (*graph_)(aLin, aCol);
}

TransitionProxy DESystem::operator()(cldes_size_t const &aLin,
                                     cldes_size_t const &aCol) {
    return TransitionProxy(this, aLin, aCol);
}

void DESystem::UpdateGraphCache_() { is_cache_outdated_ = false; }

template <typename StatesType>
DESystem::StatesSet *DESystem::Bfs_(
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
DESystem::StatesSet *DESystem::Bfs_<cldes_size_t>(
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

DESystem::StatesSet *DESystem::Bfs_() { return Bfs_(init_state_, nullptr); };

DESystem::StatesSet *DESystem::BfsCalc_(
    StatesVector &aHostX,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit,
    std::vector<cldes_size_t> const *const aStatesMap) {
    cl_uint n_initial_nodes = aHostX.size2();

    ublas::compressed_matrix<ScalarType> graph{ublas::trans(*graph_)};

    // Sum it to identity and set all nnz to 1.0f
    for (auto i = graph.begin1(); i != graph.end1(); ++i) {
        graph(i.index1(), i.index1()) = 1.0f;
    }

    // Executes BFS
    StatesVector y;
    auto n_accessed_states = 0;
    for (auto i = 0; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = ublas::prod(graph, aHostX);

        if (n_accessed_states == y.nnz()) {
            break;
        } else {
            n_accessed_states = y.nnz();
        }

        aHostX = y;
    }

    // Add results to a std::set vector
    if (aBfsVisit) {
        for (auto node = y.begin1(); node != y.end1(); ++node) {
            for (auto elem = node.begin(); elem != node.end(); ++elem) {
                aBfsVisit((*aStatesMap)[elem.index2()], node.index1());
            }
        }
        return nullptr;
    }

    auto accessed_states = new StatesSet[n_initial_nodes];
    for (auto node = y.begin1(); node != y.end1(); ++node) {
        for (auto elem = node.begin(); elem != node.end(); ++elem) {
            accessed_states[elem.index2()].emplace(node.index1());
        }
    }
    return accessed_states;
}

DESystem::StatesSet DESystem::AccessiblePart() {
    // Executes a BFS on graph_
    auto paccessible_states = Bfs_();

    auto accessible_states = *paccessible_states;
    delete[] paccessible_states;

    return accessible_states;
}

DESystem::StatesSet DESystem::CoaccessiblePart() {
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

DESystem::StatesSet DESystem::TrimStates() {
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
    // Guarantee that init_state_ is a searching node: it can be marked
    searching_nodes.emplace(init_state_);

    StatesSet accpart;
    StatesSet coaccpart = marked_states_;
    auto paccessed_states = Bfs_(
        searching_nodes,
        [this, &accpart, &coaccpart](cldes_size_t const &aInitialState,
                                     cldes_size_t const &aAccessedState) {
            if (aInitialState == this->init_state_) {
                accpart.emplace(aAccessedState);
            }

            bool const is_marked = this->marked_states_.find(aAccessedState) !=
                                   this->marked_states_.end();
            if (is_marked) {
                coaccpart.emplace(aInitialState);
            }

        });
    delete[] paccessed_states;

    StatesSet trimstates;
    std::set_intersection(accpart.begin(), accpart.end(), coaccpart.begin(),
                          coaccpart.end(),
                          std::inserter(trimstates, trimstates.begin()));

    return trimstates;
}

DESystem DESystem::Trim(bool const &aDevCacheEnabled) {
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

    DESystem trim_system{trim_graph, trimstates.size(), init_state_,
                         trim_marked_states, aDevCacheEnabled};
    return trim_system;
}

void DESystem::InsertEvents(DESystem::EventsSet const &aEvents) {
    events_ = EventsSet(aEvents);
}
