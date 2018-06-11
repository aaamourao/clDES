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
 Description: DESystem<NEvents> class implementation. DESystem<NEvents> is a
 graph, which is modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#include "cldes/TransitionProxy.hpp"
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <functional>
#include <vector>

template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>::DESystem(
  StorageIndex const& aStatesNumber,
  StorageIndex const& aInitState,
  StatesSet& aMarkedStates,
  bool const& aDevCacheEnabled)
{
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
DESystem<NEvents>::DESystem(DESystem<NEvents> const &aSys) {
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

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::GraphHostData
cldes::DESystem<NEvents, StorageIndex>::GetGraph() const
{
    return graph_;
}

template<size_t NEvents, typename StorageIndex>
void
cldes::DESystem<NEvents, StorageIndex>::CacheGraph_()
{
    is_cache_outdated_ = false;
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::EventsSet const
cldes::DESystem<NEvents, StorageIndex>::operator()(
  StorageIndex const& aLin,
  StorageIndex const& aCol) const
{
    return graph_.coeff(aLin, aCol);
}

template<size_t NEvents, typename StorageIndex>
cldes::TransitionProxy<NEvents, StorageIndex>
cldes::DESystem<NEvents, StorageIndex>::operator()(StorageIndex const& aLin,
                                                   StorageIndex const& aCol)
{
    return TransitionProxy<NEvents, StorageIndex>(this, aLin, aCol);
}

template<size_t NEvents, typename StorageIndex>
void
cldes::DESystem<NEvents, StorageIndex>::UpdateGraphCache_()
{
    is_cache_outdated_ = false;
}

template<size_t NEvents, typename StorageIndex>
template<typename StatesType>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet*
cldes::DESystem<NEvents, StorageIndex>::Bfs_(
  StatesType const& aInitialNodes,
  std::function<void(StorageIndex const&, StorageIndex const&)> const&
    aBfsVisit)
{
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    // There is no need of search if a marked state is coaccessible
    StatesVector host_x{ static_cast<StorageIndex>(states_number_),
                         static_cast<StorageIndex>(aInitialNodes.size()) };

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    std::vector<StorageIndex> states_map;
    for (auto state : aInitialNodes) {
        host_x.coeffRef(state, states_map.size()) = true;
        // Maping each search from each initial node to their correspondent
        // vector on the matrix
        states_map.push_back(state);
    }

    return BfsCalc_(host_x, aBfsVisit, &states_map);
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet*
cldes::DESystem<NEvents, StorageIndex>::Bfs_(
  StorageIndex const& aInitialNode,
  std::function<void(StorageIndex const&, StorageIndex const&)> const&
    aBfsVisit)
{
    /*
     * BFS on a Linear Algebra approach:
     *     Y = G^T * X
     */
    StatesVector host_x{ static_cast<StorageIndex>(states_number_), 1 };

    // GPUs does not allow dynamic memory allocation. So, we have
    // to set X on host first.
    host_x.coeffRef(aInitialNode, 0) = true;

    return BfsCalc_(host_x, aBfsVisit, nullptr);
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet*
cldes::DESystem<NEvents, StorageIndex>::Bfs_()
{
    return Bfs_(init_state_, nullptr);
};

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet*
cldes::DESystem<NEvents, StorageIndex>::BfsCalc_(
  StatesVector& aHostX,
  std::function<void(StorageIndex const&, StorageIndex const&)> const&
    aBfsVisit,
  std::vector<StorageIndex> const* const aStatesMap)
{
    using RowIterator =
      Eigen::InnerIterator<DESystem<NEvents, StorageIndex>::StatesVector>;

    StorageIndex n_initial_nodes = aHostX.cols();

    // Executes BFS
    StatesVector y{ static_cast<StorageIndex>(states_number_),
                    static_cast<StorageIndex>(n_initial_nodes) };
    StorageIndex n_accessed_states = 0l;
    for (StorageIndex i = 0ul; i < states_number_; ++i) {
        // Using auto bellow results in compile error
        // on the following for statement
        y = bit_graph_ * aHostX;

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        aHostX = y;
    }

    // Add results to a std::set vector
    if (aBfsVisit) {
        for (StorageIndex s = 0l; s < y.rows(); ++s) {
            for (RowIterator e(y, s); e; ++e) {
                aBfsVisit((*aStatesMap)[e.col()], e.row());
            }
        }
        return nullptr;
    }

    auto accessed_states = new StatesSet[n_initial_nodes];
    for (StorageIndex s = 0l; s < y.cols(); ++s) {
        for (RowIterator e(y, s); e; ++e) {
            accessed_states[e.col()].emplace(e.row());
        }
    }
    return accessed_states;
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet
cldes::DESystem<NEvents, StorageIndex>::AccessiblePart()
{
    // Executes a BFS on graph_
    auto paccessible_states = Bfs_();

    auto accessible_states = *paccessible_states;
    delete[] paccessible_states;

    return accessible_states;
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet
cldes::DESystem<NEvents, StorageIndex>::CoaccessiblePart()
{
    using RowIteratorConst =
      Eigen::InnerIterator<DESystem<NEvents, StorageIndex>::StatesVector const>;

    DESystem<NEvents, StorageIndex>::BitGraphHostData const invgraph =
      bit_graph_.transpose();

    StorageIndex const n_marked =
      static_cast<StorageIndex>(marked_states_.size());
    StatesVector x{ static_cast<StorageIndex>(states_number_), n_marked };
    x.reserve(marked_states_.size());

    {
        auto pos = 0ul;
        for (auto state : marked_states_) {
            x.coeffRef(state, pos) = true;
            ++pos;
        }
    }

    StatesVector y{ static_cast<StorageIndex>(states_number_), n_marked };
    auto n_accessed_states = 0l;
    for (StorageIndex i = 0; i < states_number_; ++i) {
        y = invgraph * x;

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        x = y;
    }

    y.pruned();

    StatesSet coaccessible_states;
    for (StorageIndex s = 0; s < y.outerSize(); ++s) {
        for (RowIteratorConst e(y, s); e; ++e) {
            coaccessible_states.emplace(e.row());
        }
    }

    return coaccessible_states;
}

template<size_t NEvents, typename StorageIndex>
typename cldes::DESystem<NEvents, StorageIndex>::StatesSet
cldes::DESystem<NEvents, StorageIndex>::TrimStates()
{
    using RowIteratorConst =
      Eigen::InnerIterator<DESystem<NEvents, StorageIndex>::StatesVector const>;
    using BitTriplet = Eigen::Triplet<bool, StorageIndex>;

    StatesSet const accpartstl = AccessiblePart();
    spp::sparse_hash_set<StorageIndex> accpart;
    for (StorageIndex s : accpartstl) {
        accpart.insert(s);
    }

    DESystem<NEvents, StorageIndex>::BitGraphHostData const invgraph =
      bit_graph_.transpose();

    StorageIndex const n_marked =
      static_cast<StorageIndex>(marked_states_.size());

    StatesVector x{ static_cast<StorageIndex>(states_number_), n_marked };
    x.reserve(marked_states_.size());
    std::vector<BitTriplet> xtriplet;

    {
        StorageIndex pos = 0;
        for (StorageIndex state : marked_states_) {
            xtriplet.push_back(BitTriplet(state, pos, true));
            ++pos;
        }
    }
    x.setFromTriplets(xtriplet.begin(),
                      xtriplet.end(),
                      [](bool const&, bool const&) { return true; });

    StatesVector y{ static_cast<StorageIndex>(states_number_), n_marked };
    StorageIndex n_accessed_states = 0l;
    for (StorageIndex i = 0; i < states_number_; ++i) {
        y = invgraph * x;

        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }

        x = y;
    }

    y.pruned();

    StatesSet trimstates;
    for (StorageIndex s = 0; s < y.outerSize(); ++s) {
        for (RowIteratorConst e(y, s); e; ++e) {
            if (accpart.find(e.row()) != accpart.end()) {
                trimstates.emplace(e.row());
            }
        }
    }

    return trimstates;
}

template<size_t NEvents, typename StorageIndex>
void
cldes::DESystem<NEvents, StorageIndex>::Trim()
{
    using RowIteratorGraph =
      Eigen::InnerIterator<DESystem<NEvents, StorageIndex>::GraphHostData>;
    using Triplet = Eigen::Triplet<EventsSet, StorageIndex>;
    using BitTriplet = Eigen::Triplet<bool, StorageIndex>;

    auto trimstates = this->TrimStates();

    if (trimstates.size() == static_cast<size_t>(graph_.rows())) {
        return;
    }

    // States map: old state pos -> new state pos
    std::vector<StorageIndex> statesmap(states_number_, -1);

    // Copy graph and resize it
    auto const old_graph = graph_;
    states_number_ = trimstates.size();
    graph_.resize(static_cast<StorageIndex>(states_number_),
                  static_cast<StorageIndex>(states_number_));
    bit_graph_.resize(static_cast<StorageIndex>(states_number_),
                      static_cast<StorageIndex>(states_number_));

    states_events_.erase(states_events_.begin() + states_number_,
                         states_events_.end());
    inv_states_events_.erase(inv_states_events_.begin() + states_number_,
                             inv_states_events_.end());

    events_.reset();

    // Calculate the sparsity pattern
    StorageIndex sparcitypattern = events_.count() * states_number_;

    // Fill statesmap
    {
        StorageIndex d = 0;
        for (StorageIndex s : trimstates) {
            statesmap[s] = d;
            ++d;
        }
    }

    std::vector<Triplet> triplet;
    std::vector<BitTriplet> bittriplet;

    triplet.reserve(sparcitypattern);
    bittriplet.reserve(sparcitypattern);

    // Build new graph_ slice by slice
    {
        StorageIndex row_id = 0;
        for (StorageIndex s : trimstates) {
            if (states_events_.size() > 0) {
                states_events_[row_id].reset();
                inv_states_events_[row_id].reset();
            }
            for (RowIteratorGraph e(old_graph, s); e; ++e) {
                if (statesmap[e.col()] != -1) {
                    StorageIndex const col_id = statesmap[e.col()];

                    triplet.push_back(Triplet(row_id, col_id, e.value()));
                    bittriplet.push_back(BitTriplet(col_id, row_id, true));
                    events_ |= e.value();
                    if (states_events_.size() > 0) {
                        states_events_[row_id] |= e.value();
                        inv_states_events_[col_id] |= e.value();
                    }
                }
            }
            ++row_id;
        }
    }

    // Remove aditional space
    graph_.setFromTriplets(triplet.begin(), triplet.end());
    bit_graph_.setFromTriplets(bittriplet.begin(),
                               bittriplet.end(),
                               [](bool const&, bool const&) { return true; });

    graph_.makeCompressed();
    bit_graph_.makeCompressed();

    // Calculate new marked states
    auto const old_marked = marked_states_;
    marked_states_.clear();
    for (StorageIndex s : old_marked) {
        if (statesmap[s] != -1) {
            marked_states_.emplace(statesmap[s]);
        }
    }

    return;
}

template<size_t NEvents, typename StorageIndex>
void
cldes::DESystem<NEvents, StorageIndex>::InsertEvents(
  cldes::DESystem<NEvents, StorageIndex>::EventsSet const& aEvents)
{
    events_ = EventsSet(aEvents);
}