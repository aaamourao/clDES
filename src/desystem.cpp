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
#include <functional>
#include <vector>
#include "des/transition_proxy.hpp"

using namespace cldes;

DESystem::DESystem(cldes_size_t const &aStatesNumber,
                   cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                   bool const &aDevCacheEnabled) {
    init_state_ = aInitState;
    states_number_ = aStatesNumber;
    marked_states_ = aMarkedStates;
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;

    // Resize graphs and do not preserve elements
    graph_.resize(states_number_, states_number_);
    bit_graph_.resize(states_number_, states_number_);

    // Change graphs storage type to CSR
    graph_.makeCompressed();
    bit_graph_.makeCompressed();

    // Initialize bit graph with Identity
    bit_graph_.setIdentity();

    // Reserve mem for hash tables for avoiding re-hashing
    states_events_ = StatesEventsTable(states_number_);
    inv_states_events_ = StatesEventsTable(states_number_);

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        this->CacheGraph_();
    }
}

/*
DESystem::DESystem(DESystem const &aSys) {
    init_state_ = cldes_size_t{aSys.init_state_};
    states_number_ = cldes_size_t{aSys.states_number_};
    marked_states_ = StatesSet{aSys.marked_states_};
    dev_cache_enabled_ = bool{aSys.dev_cache_enabled_};
    is_cache_outdated_ = bool{aSys.is_cache_outdated_};
    events_ = EventsSet{aSys.events_};
    graph_ = GraphHostData{aSys.graph_};
    bit_graph_ = BitGraphHostData{aSys.bit_graph_};
    states_events_ = StatesEventsTable{aSys.states_events_};
    inv_states_events_ = StatesEventsTable{aSys.inv_states_events_};
}
*/

DESystem::GraphHostData DESystem::GetGraph() const { return graph_; }

void DESystem::CacheGraph_() { is_cache_outdated_ = false; }

DESystem::EventsSet const DESystem::operator()(cldes_size_t const &aLin,
                                               cldes_size_t const &aCol) const {
    return graph_.coeff(aLin, aCol);
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
    StatesVector host_x{static_cast<Eigen::Index>(states_number_),
                        static_cast<Eigen::Index>(aInitialNodes.size())};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    std::vector<cldes_size_t> states_map;
    for (auto state : aInitialNodes) {
        host_x.coeffRef(state, states_map.size()) = true;
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
    StatesVector host_x{static_cast<Eigen::Index>(states_number_), 1};

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x.coeffRef(aInitialNode, 0) = true;

    return BfsCalc_(host_x, aBfsVisit, nullptr);
}

DESystem::StatesSet *DESystem::Bfs_() { return Bfs_(init_state_, nullptr); };

using RowIteratorConst = Eigen::InnerIterator<DESystem::StatesVector const>;
using RowIterator = Eigen::InnerIterator<DESystem::StatesVector>;

DESystem::StatesSet *DESystem::BfsCalc_(
    StatesVector &aHostX,
    std::function<void(cldes_size_t const &, cldes_size_t const &)> const
        &aBfsVisit,
    std::vector<cldes_size_t> const *const aStatesMap) {
    cldes_size_t n_initial_nodes = aHostX.cols();

    // Executes BFS
    StatesVector y{static_cast<Eigen::Index>(states_number_),
                   static_cast<Eigen::Index>(n_initial_nodes)};
    auto n_accessed_states = 0l;
    for (auto i = 0ul; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = (bit_graph_ * aHostX).pruned();

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        aHostX = y;
    }

    // Add results to a std::set vector
    if (aBfsVisit) {
        for (auto s = 0l; s < y.rows(); ++s) {
            for (RowIterator e(y, s); e; ++e) {
                aBfsVisit((*aStatesMap)[e.col()], e.row());
            }
        }
        return nullptr;
    }

    auto accessed_states = new StatesSet[n_initial_nodes];
    for (auto s = 0l; s < y.cols(); ++s) {
        for (RowIterator e(y, s); e; ++e) {
            accessed_states[e.col()].emplace(e.row());
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
    DESystem::BitGraphHostData const invgraph = bit_graph_.transpose();

    Eigen::Index const n_marked =
        static_cast<Eigen::Index>(marked_states_.size());
    StatesVector x{static_cast<Eigen::Index>(states_number_), n_marked};
    x.reserve(marked_states_.size());

    {
        auto pos = 0ul;
        for (auto state : marked_states_) {
            x.coeffRef(state, pos) = true;
            ++pos;
        }
    }

    StatesVector y{static_cast<Eigen::Index>(states_number_), n_marked};
    auto n_accessed_states = 0l;
    for (auto i = 0ul; i < states_number_; ++i) {
        y = (invgraph * x).pruned();

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        x = y;
    }

    y.pruned();

    StatesSet coaccessible_states;
    for (auto s = 0; s < y.outerSize(); ++s) {
        for (RowIteratorConst e(y, s); e; ++e) {
            coaccessible_states.emplace(e.row());
        }
    }

    return coaccessible_states;
}

DESystem::StatesSet DESystem::TrimStates() {
    StatesSet const accpartstl = AccessiblePart();
    QSet<cldes_size_t> accpart;
    for (auto s : accpartstl) {
        accpart.insert(s);
    }

    DESystem::BitGraphHostData const invgraph = bit_graph_.transpose();

    Eigen::Index const n_marked =
        static_cast<Eigen::Index>(marked_states_.size());
    StatesVector x{static_cast<Eigen::Index>(states_number_), n_marked};
    x.reserve(marked_states_.size());

    {
        auto pos = 0ul;
        for (auto state : marked_states_) {
            x.coeffRef(state, pos) = true;
            ++pos;
        }
    }

    StatesVector y{static_cast<Eigen::Index>(states_number_), n_marked};
    auto n_accessed_states = 0l;
    for (auto i = 0ul; i < states_number_; ++i) {
        y = (invgraph * x).pruned();

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        x = y;
    }

    y.pruned();

    StatesSet trimstates;
    for (auto s = 0; s < y.outerSize(); ++s) {
        for (RowIteratorConst e(y, s); e; ++e) {
            if (accpart.contains(e.row())) {
                trimstates.emplace(e.row());
            }
        }
    }

    return trimstates;
}

using RowIteratorGraph = Eigen::InnerIterator<DESystem::GraphHostData>;

void DESystem::Trim() {
    auto trimstates = this->TrimStates();

    if (trimstates.size() == static_cast<unsigned long>(graph_.rows())) {
        return;
    }

    // Copy graph and resize it
    auto old_graph = graph_;
    auto const old_states_number = states_number_;
    states_number_ = trimstates.size();
    graph_.resize(static_cast<long>(states_number_),
                  static_cast<long>(states_number_));
    bit_graph_.resize(static_cast<long>(states_number_),
                      static_cast<long>(states_number_));

    // States map: old state pos -> new state pos
    std::vector<long> statesmap(old_states_number, -1);

    states_events_.erase(states_events_.begin() + states_number_,
                         states_events_.end());
    inv_states_events_.erase(inv_states_events_.begin() + states_number_,
                             inv_states_events_.end());

    // Calculate the sparsity pattern
    auto sparcitypattern = 0ul;
    for (auto sit = trimstates.begin(); sit != trimstates.end(); ++sit) {
        auto const d = std::distance(trimstates.begin(), sit);
        sparcitypattern += old_graph.row(*sit).nonZeros();
        statesmap[*sit] = d;
    }

    using Triplet = Eigen::Triplet<EventsBitArray>;
    using BitTriplet = Eigen::Triplet<bool>;

    std::vector<Triplet> triplet;
    std::vector<BitTriplet> bittriplet;

    triplet.reserve(sparcitypattern);
    bittriplet.reserve(sparcitypattern);

    // Build new graph_ slice by slice
    for (auto st = trimstates.begin(); st != trimstates.end(); ++st) {
        auto const row_id = std::distance(trimstates.begin(), st);

        for (RowIteratorGraph e(old_graph, *st); e; ++e) {
            if (statesmap[e.col()] != -1) {
                auto const col_id = statesmap[e.col()];

                triplet.push_back(Triplet(row_id, col_id, e.value()));
                bittriplet.push_back(BitTriplet(col_id, row_id, true));
                events_ |= e.value();
                states_events_[row_id] |= e.value();
                inv_states_events_[col_id] |= e.value();
            }
        }
    }

    // Remove aditional space
    graph_.setFromTriplets(triplet.begin(), triplet.end());
    bit_graph_.setFromTriplets(bittriplet.begin(), bittriplet.end());

    // Calculate new marked states
    auto const old_marked = marked_states_;
    marked_states_.clear();
    for (auto s : old_marked) {
        if (statesmap[s] != -1) {
            marked_states_.emplace(statesmap[s]);
        }
    }

    return;
}

void DESystem::InsertEvents(DESystem::EventsSet const &aEvents) {
    events_ = EventsSet(aEvents);
}
