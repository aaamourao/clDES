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

 File: cldes/src/des/DESystemCore.hpp
 description: DESystem class methods definitions
 =========================================================================
*/
/*!
 * \file cldes/src/des/DESystemCore.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-06-16
 *
 * DESystem template class definition.
 */

namespace cldes {
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>::DESystem()
{
    inv_graph_ = nullptr;
    is_cache_outdated_ = false;
    this->marked_states_ = StatesSet{};
    graph_ = GraphHostData{};
}

template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>::DESystem(StorageIndex const& aStatesNumber,
                                          StorageIndex const& aInitState,
                                          StatesSet& aMarkedStates,
                                          bool const& aDevCacheEnabled)
  : DESystemBase{ aStatesNumber, aInitState }
{
    dev_cache_enabled_ = aDevCacheEnabled;
    is_cache_outdated_ = true;
    inv_graph_ = nullptr;

    this->marked_states_ = aMarkedStates;

    // Resize graphs and do not preserve elements
    graph_ = GraphHostData{ static_cast<StorageIndexSigned>(aStatesNumber),
                            static_cast<StorageIndexSigned>(aStatesNumber) };

    // Change graphs storage type to CSR
    graph_.makeCompressed();

    // Reserve memory to make insertions efficient
    this->states_events_ = StatesEventsTable(aStatesNumber);
    this->inv_states_events_ = StatesEventsTable(aStatesNumber);

    // If device cache is enabled, cache it
    if (dev_cache_enabled_) {
        this->cacheGraph_();
    }
}

template<uint8_t NEvents, typename StorageIndex>
TransitionProxy<NEvents, StorageIndex>
DESystem<NEvents, StorageIndex>::operator()(StorageIndex const& aQfrom,
                                            StorageIndex const& aQto)
{
    return TransitionProxy<NEvents, StorageIndex>(*this, aQfrom, aQto);
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystem<NEvents, StorageIndex>::StatesSet
DESystem<NEvents, StorageIndex>::accessiblePart() const noexcept
{
    // Executes a BFS on graph_
    auto accessible_states = bfs_();
    return *accessible_states;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystem<NEvents, StorageIndex>::StatesSet
DESystem<NEvents, StorageIndex>::coaccessiblePart() const noexcept
{
    GraphHostData searchgraph{ this->states_number_, this->states_number_ };
    searchgraph.setIdentity();
    searchgraph += graph_;
    StatesVector const invgraph{ searchgraph.template cast<bool>() };

    StatesVector x{ static_cast<StorageIndexSigned>(this->states_number_), 1 };
    x.reserve(this->marked_states_.size());
    for (auto state : this->marked_states_) {
        x.coeffRef(state, 0) = true;
    }

    return *procStVec_(bfsCalc_(std::move(x), std::move(invgraph)));
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystem<NEvents, StorageIndex>::StatesSet
DESystem<NEvents, StorageIndex>::trimStates() const noexcept
{
    auto accpartstl = accessiblePart();
    spp::sparse_hash_set<StorageIndex> accpart;
    for (StorageIndex s : accpartstl) {
        accpart.insert(s);
    }
    GraphHostData ident{ static_cast<Eigen::Index>(this->states_number_),
                         static_cast<Eigen::Index>(this->states_number_) };
    ident.setIdentity();
    StatesVector const searchgraph{ (graph_ + ident).template cast<bool>() };

    StatesVector x{ static_cast<StorageIndexSigned>(this->states_number_), 1 };
    x.reserve(this->marked_states_.size());
    std::vector<BitTriplet> xtriplet;
    for (StorageIndex state : this->marked_states_) {
        xtriplet.push_back(BitTriplet(state, 0, true));
    }
    x.setFromTriplets(xtriplet.begin(),
                      xtriplet.end(),
                      [](bool const&, bool const&) { return true; });

    // StatesVector y = bfsCalc_(x, searchgraph);
    return procStVec_(bfsCalc_(std::move(x), std::move(searchgraph)),
                      [atable = std::move(accpart)](
                        StorageIndex const&, StorageIndex const& qto) -> bool {
                          if (atable.find(qto) != atable.end()) {
                              return true;
                          }
                          return false;
                      });
}

template<uint8_t NEvents, typename StorageIndex>
DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>>&
DESystem<NEvents, StorageIndex>::trim() noexcept
{
    auto trimstates = this->trimStates();
    if (trimstates.size() == static_cast<size_t>(graph_.rows())) {
        return *this;
    }
    // States map: old state pos -> new state pos
    std::vector<StorageIndexSigned> statesmap(this->states_number_, -1);
    // Copy graph and resize it
    auto const old_graph = graph_;
    this->states_number_ = trimstates.size();
    graph_.resize(static_cast<StorageIndexSigned>(this->states_number_),
                  static_cast<StorageIndexSigned>(this->states_number_));
    if (this->states_events_.size() > 0) {
        this->states_events_.erase(this->states_events_.begin() +
                                     this->states_number_,
                                   this->states_events_.end());
        this->inv_states_events_.erase(this->inv_states_events_.begin() +
                                         this->states_number_,
                                       this->inv_states_events_.end());
    }
    this->events_.reset();
    // Calculate the sparsity pattern
    // auto sparcitypattern = this->events_.count() * this->states_number_;
    unsigned long sparcitypattern = 0;
    // Fill statesmap
    {
        auto d = 0ul;
        for (StorageIndex s : trimstates) {
            statesmap[s] = d;
            ++d;
            sparcitypattern += old_graph.innerVector(s).nonZeros();
        }
    }
    std::vector<Triplet> triplet;
    triplet.reserve(sparcitypattern);
    // Build new graph_ slice by slice
    {
        auto row_id = 0ul;
        for (StorageIndex s : trimstates) {
            if (this->states_events_.size() > 0) {
                this->states_events_[row_id].reset();
                this->inv_states_events_[row_id].reset();
            }
            for (RowIterator e(old_graph, s); e; ++e) {
                if (statesmap[e.col()] != -1) {
                    auto const col_id = statesmap[e.col()];
                    triplet.push_back(Triplet(row_id, col_id, e.value()));
                    this->events_ |= e.value();
                    if (this->states_events_.size() > 0) {
                        this->states_events_[row_id] |= e.value();
                        this->inv_states_events_[col_id] |= e.value();
                    }
                }
            }
            ++row_id;
        }
    }
    // Remove aditional space
    graph_.setFromTriplets(triplet.begin(), triplet.end());
    graph_.makeCompressed();
    // Calculate new marked states
    auto const old_marked = this->marked_states_;
    this->marked_states_.clear();
    for (StorageIndex s : old_marked) {
        if (statesmap[s] != -1) {
            this->marked_states_.emplace(statesmap[s]);
        }
    }
    return *this;
}

template<uint8_t NEvents, typename StorageIndex>
inline void
DESystem<NEvents, StorageIndex>::cropGraph_() noexcept
{
    if (true == true) {
    }
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystem<NEvents, StorageIndex>::insertEvents(EventsSet const& aEvents) noexcept
{
    this->events_ = EventsSet(aEvents);
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystem<NEvents, StorageIndex>::StorageIndexSigned
DESystem<NEvents, StorageIndex>::trans_impl(StorageIndex const& aQ,
                                            ScalarType const& aEvent) const
  noexcept
{
    if (!this->states_events_[aQ].test(aEvent)) {
        return -1;
    }

    for (RowIterator qiter(graph_, aQ); qiter; ++qiter) {
        if (qiter.value().test(aEvent)) {
            return qiter.col();
        }
    }

    // It will never reach here, but the warning is bothering me
    return -1;
}

template<uint8_t NEvents, typename StorageIndex>
StatesArray<StorageIndex>
DESystem<NEvents, StorageIndex>::invtrans_impl(StorageIndex const& aQ,
                                               ScalarType const& aEvent) const
{
    StatesArray<StorageIndex> inv_trans;

    if (!this->inv_states_events_[aQ].test(aEvent)) {
        return inv_trans;
    }

    for (RowIterator qiter(*inv_graph_, aQ); qiter; ++qiter) {
        if (qiter.value().test(aEvent)) {
            inv_trans.push_back(qiter.col());
        }
    }

    // It will never reach here, but the warning is bothering me
    return inv_trans;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystem<NEvents, StorageIndex>::allocateInvertedGraph_impl() const noexcept
{
    inv_graph_ = std::make_shared<GraphHostData>(graph_.transpose());
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystem<NEvents, StorageIndex>::clearInvertedGraph_impl() const noexcept
{
    inv_graph_ = nullptr;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystem<NEvents, StorageIndex>::cacheGraph_() noexcept
{
    is_cache_outdated_ = false;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystem<NEvents, StorageIndex>::updateGraphCache_() noexcept
{
    is_cache_outdated_ = false;
}

template<uint8_t NEvents, typename StorageIndex>
template<class StatesType>
inline std::shared_ptr<typename DESystem<NEvents, StorageIndex>::StatesSet>
DESystem<NEvents, StorageIndex>::bfs_(StatesType const&& aInitialNodes) const
  noexcept
{
    // There is no need of search if a marked state is coaccessible
    StatesVector host_x{ static_cast<StorageIndexSigned>(this->states_number_),
                         static_cast<StorageIndexSigned>(
                           aInitialNodes.size()) };
    std::vector<StorageIndex> states_map;
    for (auto state : aInitialNodes) {
        host_x.coeffRef(state, states_map.size()) = true;
        // Maping each search from each initial node to their correspondent
        // vector on the matrix
        states_map.push_back(state);
    }
    GraphHostData searchgraph{ static_cast<long>(this->states_number_),
                               static_cast<long>(this->states_number_) };
    searchgraph.setIdentity();
    searchgraph += graph_;
    StatesVector const invgraph{
        searchgraph.template cast<bool>().transpose()
    };
    return procStVec_(bfsCalc_(std::move(host_x), std::move(invgraph)));
}

template<uint8_t NEvents, typename StorageIndex>
inline std::shared_ptr<typename DESystem<NEvents, StorageIndex>::StatesSet>
DESystem<NEvents, StorageIndex>::bfs_() const noexcept
{
    StatesVector host_x{ static_cast<StorageIndexSigned>(this->states_number_),
                         1 };
    host_x.coeffRef(this->init_state_, 0) = true;
    GraphHostData searchgraph{ static_cast<long>(this->states_number_),
                               static_cast<long>(this->states_number_) };
    searchgraph.setIdentity();
    searchgraph += graph_;
    StatesVector const invgraph{
        searchgraph.template cast<bool>().transpose()
    };
    return procStVec_(bfsCalc_(std::move(host_x), std::move(invgraph)));
}

template<uint8_t NEvents, typename StorageIndex>
template<class GraphT>
inline typename DESystem<NEvents, StorageIndex>::StatesVector
DESystem<NEvents, StorageIndex>::bfsCalc_(StatesVector&& aHostX,
                                          GraphT const&& aSearchGraph) const
  noexcept
{
    /*!
     * BFS on a Linear Algebra approach:
     *     \f$Y = G^T * X\f$
     */
    StatesVector y{ static_cast<StorageIndexSigned>(this->states_number_),
                    static_cast<StorageIndexSigned>(aHostX.cols()) };
    auto n_accessed_states = 0l;
    for (StorageIndex i = 0ul; i < this->states_number_; ++i) {
        y = aSearchGraph * aHostX;
        if (n_accessed_states == y.nonZeros()) {
            break;
        } else {
            n_accessed_states = y.nonZeros();
        }
        aHostX = y;
    }
    y.pruned();
    return y;
}

template<uint8_t NEvents, typename StorageIndex>
inline std::shared_ptr<typename DESystem<NEvents, StorageIndex>::StatesSet>
DESystem<NEvents, StorageIndex>::procStVec_(StatesVector const&& aY) const
  noexcept
{
    // Unfortunatelly, only C++17 allows shared_ptr to arrays
    std::shared_ptr<StatesSet> accessed_states{
        new StatesSet[aY.cols()], std::default_delete<StatesSet[]>()
    };
    for (auto s = 0l; s < aY.outerSize(); ++s) {
        for (BitIteratorConst e(aY, s); e; ++e) {
            accessed_states.get()[e.col()].emplace(e.row());
        }
    }
    return accessed_states;
}

template<uint8_t NEvents, typename StorageIndex>
template<typename FuncT>
inline typename DESystem<NEvents, StorageIndex>::StatesSet
DESystem<NEvents, StorageIndex>::procStVec_(
  StatesVector const&& aY,
  FuncT const&& aF,
  std::shared_ptr<std::vector<StorageIndex> const> const& aStatesMap) const
  noexcept
{
    StatesSet processedstates{};
    if (aStatesMap) {
        for (auto s = 0l; s < aY.outerSize(); ++s) {
            for (BitIteratorConst e(aY, s); e; ++e) {
                if (aF((*aStatesMap)[e.col()], e.row())) {
                    processedstates.emplace(e.row());
                }
            }
        }
    } else {
        for (auto s = 0l; s < aY.outerSize(); ++s) {
            for (BitIteratorConst e(aY, s); e; ++e) {
                if (aF(e.col(), e.row())) {
                    processedstates.emplace(e.row());
                }
            }
        }
    }
    return processedstates;
}

template<uint8_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>&
DESystem<NEvents, StorageIndex>::proj_impl(EventsSet const& aAlphabet) noexcept
{
    for (StorageIndex q = 0; q < graph_.rows(); ++q) {
        for (RowIteratorGraph d(graph_, q); d; ++d) {
            d.valueRef() &= aAlphabet;
        }
        if (this->states_events_.size() > 0) {
            this->states_events_[q] &= aAlphabet;
            this->inv_states_events_[q] &= aAlphabet;
        }
    }
    graph_.prune(EventsSet{ 0 });
    return *this;
}

template<uint8_t NEvents, typename StorageIndex>
bool
cldes::DESystem<NEvents, StorageIndex>::checkObsProp_impl(
  EventsSet const&) const noexcept
{
    return false;
}
}
