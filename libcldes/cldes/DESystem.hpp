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

 File: cldes/DESystem.hpp
 description: DESystem template class declaration and definition . DESystem is a
 graph, which is modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/
/*!
 * \file cldes/DESystem.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-06-16
 *
 * DESystem template class declaration and definition . DESystem is a
 * graph, which is modeled as a Sparce Adjacency Matrix.
 */

#ifndef DESYSTEM_HPP
#define DESYSTEM_HPP

#include "cldes/DESystemBase.hpp"
#include "cldes/src/des/DESystemFwd.hpp"

namespace cldes {

/*! \class DESystem
 *  \brief A discrete-events system on host memory
 *  \details DESystem implements a discrete-event system as a graph stored on
 *  a sparse adjacency matrix.
 *
 *  Single systems operations are available as member functions:
 *
 *  |    Operation         |       Method           |
 *  |:--------------------:|:----------------------:|
 *  | Accessible part      | accessiblePart()       |
 *  | Coaccessible part    | coaccessiblePart()     |
 *  | Get trim states      | trimStates()           |
 *  | trim                 | trim()                 |
 *  | transition           | trans()                |
 *  | invtransition        | invtrans()             |
 *
 *  Where S is the number of states, E the number of events and T the number
 *  of transitions (considering a events OR expression a single transition)
 *
 * \tparam NEvents Number of events: max 255
 * \tparam StorageIndex Unsigned type used for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystem
  : public DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>>
{
public:
    /*! \brief Base alias
     * \details Alias to base class with implicit template params
     */
    using DESystemBase =
      DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>>;

    /*! \brief StorageIndex signed type
     * \details Eigen uses signed indexes
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief EventsSet
     *  \details Set containing 8bit intergets which represent events.
     */
    using EventsSet = EventsSet<NEvents>;

    /*! \brief Adjacency matrix of bitarrays implementing a directed graph
     * \details The graph represents the DES automata:
     * Non zero element: contains at least one transition
     * Each non 0 bit of each element: event that lead to the next
     * stage
     * * row index: from state
     * * col index: to state
     */
    using GraphHostData = Eigen::SparseMatrix<EventsSet, Eigen::RowMajor>;

    /*! \brief Set of states type
     *  \details Set containg unsigned interget types which represent states.
     */
    using StatesSet = typename DESystemBase::StatesSet;

    /*! \brief Table of transitions on a STL container
     *  \details Vector containing bitsets which represent the events that a
     *  state represented by the vector index contains.
     */
    using StatesEventsTable = typename DESystemBase::StatesEventsTable;

    /*! \brief Bit adjacency matrix implementing a directed graph
     * \details Sparse matrix implementing a graph with an adjacency matrix.
     * The graph represents a simplified DES automata:
     * * Non zero element: contains at least one transition
     * * row index: from state
     * * col index: to state
     */
    using BitGraphHostData = Eigen::SparseMatrix<bool, Eigen::ColMajor>;
    using BitRowGraphHostData = Eigen::SparseMatrix<bool, Eigen::RowMajor>;

    /*! \brief Adjacency matrix of bit implementing searching nodes
     * \details Structure used for traversing the graph using a linear algebra
     * approach
     */
    using StatesDenseVector = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
    using StatesVector = Eigen::SparseMatrix<bool, Eigen::ColMajor>;

    /*! \brief Set of Events implemented as a Hash Table for searching
     * efficiently.
     *  \details Sparse hash set used to encapsulate events (8bit intergers)
     */
    using EventsTable = spp::sparse_hash_set<uint8_t>;

    /*! \brief Vector of states type
     *  \details Vector of usigned interger type which represent states.
     */
    using StatesTable = typename DESystemBase::StatesTable;

    /*! \brief Vector of inverted transitions
     * \ details f(s, e) = s_out -> this vector stores = (s_out, (s, e))
     */
    using TrVector =
      std::vector<std::pair<StorageIndex, InvArgtrans<StorageIndex>*>>;

    /*! \brief Graph const iterator
     * \details Used to iterate over the bfs result
     */
    using ColIteratorConst = Eigen::InnerIterator<StatesVector const>;
    using DenseColIteratorConst = Eigen::InnerIterator<StatesDenseVector const>;
    using BitRowIteratorConst = Eigen::InnerIterator<BitRowGraphHostData const>;
    using BitIteratorConst = Eigen::InnerIterator<BitGraphHostData const>;

    /*! \brief Graph const iterator
     * \details Used to iterate over the adjacency matrix
     */
    using RowIterator = Eigen::InnerIterator<GraphHostData const>;

    /*! \brief Graph iterator
     * \details Used to iterate over the adjacency matrix
     */
    using RowIteratorGraph = typename GraphHostData::InnerIterator;

    /*! \brief Triplet type
     * \details 3-tuples (qfrom, qto, event) for filling adjacency matrix
     * efficiently
     */
    using Triplet = Triplet<NEvents>;

    /*! \brief Default default constructor
     *  \details Creates an empty system:
     *  * Number of states = 0
     *  * Initial state = 0
     *  * Marked states = { }
     *  * Bit graph = Identity
     *  * Graph = Empty matrix with 0 rows and 0 columns
     */
    DESystem();

    /*! \brief DESystem constructor
     * \details Create a system with certain parameters.
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     * @param aMarkedStates System's marked states
     * @param aDevCacheEnabled Enable or disable device cache for graph data
     */
    DESystem(StorageIndex const& aStatesNumber,
             StorageIndex const& aInitState,
             StatesSet& aMarkedStates,
             bool const& aDevCacheEnabled = true);

    /*! \brief DESystem destructor
     *  \details Override base destructor for avoiding memory leaks.
     */
    ~DESystem() = default;

    /*! \brief Move constructor
     *  \details Enable move semantics
     */
    DESystem(DESystem&&) = default;

    /*! \brief Copy constructor
     *  \details Enable copy by reference of a const system.
     */
    DESystem(DESystem const&) = default;

    /*! \brief Operator =
     *  \details Use move semantics with operator = when assigning to rvalues.
     */
    DESystem<NEvents, StorageIndex>& operator=(DESystem&&) = default;

    /*! \brief Operator = to const type
     *  \details Enable copy by reference of const systems with operator =.
     */
    DESystem<NEvents, StorageIndex>& operator=(DESystem const&) = default;

    /*! \brief clone method for polymorphic copy
     *  \return Shared pointer to this object
     */
    std::shared_ptr<DESystemBase> constexpr clone_impl() const noexcept
    {
        std::shared_ptr<DESystemBase> this_ptr =
          std::make_shared<DESystem>(*this);
        return this_ptr;
    }

    /*! \brief Check if this system is a virtual proxy
     *
     *  \details Operations between two systems return virtual proxies to
     *  allow lazy evaluation of large systems while it is not necessary
     *  to obtain the whole system. While the system is virtual, the
     *  memory complexity of DESystem is cheaper, and the time complexities
     *  of its embedded operations are more expensive.
     *
     *  \return True: DESystem is always a real object
     */
    bool constexpr static isVirtual_impl() noexcept { return false; }

    /*! \brief Graph getter
     *
     *  \return Eigen sparse matrix of bitset representing the sysmte on
     *  compressed mode.
     */
    GraphHostData constexpr getGraph() const noexcept { return graph_; }

    /*! \brief Returns events that lead a transition between two states
     *
     * @param aQfrom State where the transition departs.
     * @param aQto State where the transition arrives.
     * \return A bitset where each set bit index is a event that lead to the
     * next state.
     */
    EventsSet const constexpr operator()(StorageIndex const& aQfrom,
                                         StorageIndex const& aQto) const
      noexcept
    {
        return graph_.coeff(aQfrom, aQto);
    }

    /*! \brief Returns reference events that lead a transition between two
     * states
     *  \details This operator can be used for getting events from a transition
     *  between two systems or to assign it. It was implemented to keep track
     *  when the copy of its system on the device memory is outdated.
     *
     * @param aQfrom State where the transition departs.
     * @param aQto State where the transition arrives.
     * \return A reference to the event's transition.
     */
    TransitionProxy<NEvents, StorageIndex> operator()(
      StorageIndex const& aQfrom,
      StorageIndex const& aQto);

    /*! \brief Returns state set containing the accessible part of automa
     * \details Executes a breadth-first search in bit the graph, which
     * represents the DES, starting from its initial state. It returns a set
     * containing all nodes which are accessible from the initial state. This
     * operation is implemented on a linear algebra approach by a multiplication
     * between two sparse matrices.
     *
     * \return A set of the accessed states
     */
    StatesSet accessiblePart() const noexcept;

    /*! \brief Returns state set containing the coaccessible part of automata
     * \details Inverts the bit graph by transposing it (does not store it
     * on the bit_graph_ member). Makes a breadth-first search starting from
     * the marked states all at once by multiplying it by a sparse matrix with
     * states_number_ rows and marked_states_.size() columns.
     *
     * \return A set containing the accessed states by the bfs (coacessible
     * states)
     */
    StatesSet coaccessiblePart() const noexcept;

    /*! \brief Get states set which represent the trim states of the current
     * system
     * \details Gets the intersection between the accessible part and the
     * coaccessible part.
     *
     * \return A set containing the trim states intergers.
     */
    StatesSet trimStates() const noexcept;

    /*! \brief Returns DES which is the trim part of this
     * \details Cut the non-trim states from the current system.
     * TODO: Make it lazy. Maybe setting values by 0 and making the graph_
     * compressed again would be a good solution.
     *
     * \return Reference to this system
     */
    DESystemBase& trim() noexcept;

    /*! \brief Insert events
     * \details Set the member events_ with a set containing all events
     * that are present
     * on the current system.
     * \warning Take care using this method. It was designed for
     * testing and debugging.
     *
     * @param aEvents Set containing all the new events of the current system
     * \return void
     */
    void insertEvents(EventsSet const& aEvents) noexcept;

    /*! \brief Check if transition exists
     * \details given a transition f(q, e) = qto, it
     * checks if qto is different of empty.
     *
     * @param aQ State
     * @param aEvent Event
     * \return Returns true if DES transition exists, false otherwise
     */
    bool constexpr containstrans_impl(StorageIndex const& aQ,
                                      ScalarType const& aEvent) const noexcept
    {
        return this->states_events_[aQ].test(aEvent);
    }

    /*! \brief transition function
     * \details given a transition f(q, e) = qto, it
     * calculates if qto is different of empty.
     *
     * @param aQ State
     * @param aEvent Event
     * \return The state where the transition leads or -1 when it is empty
     */
    StorageIndexSigned trans_impl(StorageIndex const& aQ,
                                  ScalarType const& aEvent) const noexcept;

    /*! \brief Check if the current system contains at least one inverse
     * transition
     * \warning This method requires to run AllocateInvGraph previously. It
     * can be very expensive for large systems.
     *
     * @param aQ State
     * @param aEvent Event
     * \return True if DES inverse transition exists, false otherwise
     */
    bool constexpr containsinvtrans_impl(StorageIndex const& aQ,
                                         ScalarType const& aEvent) const
    {
        return this->inv_states_events_[aQ].test(aEvent);
    }

    /*! \brief DES Inverse transition function
     * \warning This method requires to run AllocateInvGraph previously. It
     * can be very expensive for large systems, because it implies in
     * transposing the adjacency matrix.
     *
     * @param aQ State
     * @param aEvent Event
     * \return DES states array containing the inverse transition.
     * It may be empty, if there is none.
     */
    StatesArray<StorageIndex> invtrans_impl(StorageIndex const& aQ,
                                            ScalarType const& aEvent) const;

    /*! \brief Get events of all transitions of a specific state
     * \details Since this is information is stored on a vector on concrete
     * systems, this operation is really cheap, O(1).
     *
     * @param aQ A state on the sys
     * \return Returns EventsSet relative to state q
     */
    EventsSet constexpr getStateEvents_impl(StorageIndex const& aQ) const
      noexcept
    {
        return this->states_events_[aQ];
    }

    /*! \brief Get events of all transitions that lead to a specific state
     * \details Since this is information is stored on a vector on concrete
     * systems, this operation is really cheap, O(1).
     *
     * @param aQ A state on the sys
     * \return Returns EventsSet relative to inverse transitions of the state q
     */
    EventsSet constexpr getInvStateEvents_impl(StorageIndex const& aQ) const
    {
        return this->inv_states_events_[aQ];
    }

    /*! \brief Invert graph
     * \details This is used in the monolithic supervisor synthesis. It needs
     * to be executed before making inverse transitions operations.
     * It is const, since it changes only a mutable member
     * \warning It can be very inneficient for very large graphs
     *
     * \return void
     */
    void allocateInvertedGraph_impl() const noexcept;

    /*! \brief Free inverted graph
     * It is const, since it changes only a mutable member
     * \warning It is recommended to clear the it make the inverted graph clear
     * as much as possible when handling large systems.
     *
     * \return void
     */
    void clearInvertedGraph_impl() const noexcept;

    /*! \brief Observer property checker
     *
     * @param aLang Language
     */
    bool checkObsProp_impl(EventsSet const& aAlphabet) const noexcept;

    /*! \brief Projection operation
     *
     * @param aAlphanet Events set which the system will be projected.
     */
    DESystem& proj_impl(EventsSet const& aAlphabet) noexcept;

    /*! \brief Inverse projection operation
     *
     * @param aLang Set of events
     */
    DESystem constexpr& invproj_impl(EventsSet const&) const noexcept
    {
        return *this;
    }

protected:
    /*! \brief Method for caching the graph
     * \details Copy the graph after transposing it to the device memory.
     *
     * \return void
     */
    void cacheGraph_() noexcept;

    /*! \brief Method for updating the graph
     * \details Refresh the existent graph data on device memory.
     *
     * \return void
     */
    void updateGraphCache_() noexcept;

    /*! \brief Setup BFS and return accessed states array
     * \details Executes a breadth-first search on the graph starting from N
     * nodes in aInitialNodes. The algorithm is based on SpGEMM.
     *
     * @param aInitialNodes Set of nodes where the searches will start
     * @param aBfsVisit Function to execute on every accessed state: it can be
     * null
     * \return Return a set containing the accessed states, or nullptr if none
     * was accessed or a bfs visit function was defined
     */
    template<class StatesType>
    inline std::shared_ptr<StatesSet> bfs_(
      StatesType const&& aInitialNodes) const noexcept;

    /*! \brief Overload bfs_ for the special case of starting from the inital
     * state
     * \details Executes a breadth-first search on the graph starting from the
     * system's initial node The algorithm is based on SpGEMM.
     *
     * \return Return a set containing the accessed states, or nullptr if none
     * was accessed or a bfs visit function was defined
     */
    inline std::shared_ptr<StatesSet> bfs_() const noexcept;

    /*! \brief Calculates Bfs and returns accessed states array
     * \details Executes a breadth first search on the graph starting from one
     * single node. The algorithm is based on SpGEMM.
     *
     * @param aHostX Bit vector with initial nodes as nonzeros rows.
     * @param aSearchGraph graph where the search are going to be executed
     */
    template<class GraphT>
    inline StatesVector bfsCalc_(StatesVector&& aHostX,
                                 GraphT const&& aSearchGraph) const noexcept;

    /* \brief Process result of vectorial operations
     *
     * @param aBfsVisit Function to execute on every accessed state: it can be
     * null
     * @param aStatesMap Map each vector used to implement in a bfs when
     * multiple ones are execute togheter: it can be null
     */
    inline std::shared_ptr<StatesSet> procStVec_(StatesVector const&& y) const
      noexcept;

    /* \brief Overload
     *
     */
    template<typename FuncT = std::function<bool(StorageIndex const&,
                                                 StorageIndex const&) const>>
    inline StatesSet procStVec_(
      StatesVector const&& y,
      FuncT const&& aF,
      std::shared_ptr<std::vector<StorageIndex> const> const& aStatesMap =
        nullptr) const noexcept;

    /*! \brief Remove vertices from graph_
     */
    inline void cropGraph_(
      unsigned long const& aSparcityaPattern,
      GraphHostData const&& aOldGraph,
      StatesSet const&& aTrimStates,
      std::vector<StorageIndexSigned> const& aStatesMap) noexcept;

    /*! \brief Remove marked states that are blocking
     */
    inline void cropMarkedStates_(
      std::vector<StorageIndexSigned> const&& aStatesMap) noexcept;

private:
    friend class TransitionProxy<NEvents, StorageIndex>;
#ifdef CLDES_OPENCL_ENABLED
    friend class DESystemCL<NEvents, StorageIndex>;
#endif

    /*! \brief DESystem operations virtual proxies
     *
     * Functions which implement DES operations between two systems and need
     * to access private members of DESystem. The design constrains that lead to
     * declare the following friend(s):
     *
     * * Lazy evaluation
     * * Functions which have more than one system as input and returns a
     * different system should not be a member of DESystem class.
     * * More efficiency when accessing vars than using getters and setters.
     * * Define a different namespace for DES operations.
     */
    friend class op::SyncSysProxy<NEvents, StorageIndex>;

    /*! \brief Graph represented by an adjacency matrix of bitsets
     * \details A sparse matrix which represents the automata as a graph in an
     * adjascency matrix. It is implemented as a CSR scheme.
     *  * Non zero element: events bit from <row index> to <col index> that
     * lead to starting node (row index) to the arriving node (col index).
     * * e.g. M(2, 3) = 101; transition from state 2 to state 3 with the
     * condition event 0 OR event 2.
     */
    GraphHostData graph_;

    /*! \brief Inverted graph shared pointer
     * \details Used for searching inverted transitions when necessary.
     * Keep it free whenever it is possible.
     * TODO: Make unique
     */
    std::shared_ptr<GraphHostData> mutable inv_graph_;

    /*! \brief Keeps if caching graph data on device is enabled
     * \details If dev_cache_enabled_ is true, the graph should be cached on
     * the device memory, so device_graph_ is not nullptr. It can be set at any
     * time at run time, so it is not a constant.
     * \warning The OpenCL functions are disabled.
     */
    bool dev_cache_enabled_;

    /*! \brief Keeps track if the device graph cache is outdated
     * \details  Tracks if cache needs to be updated or not.
     * \warning The OpenCL functions are disabled.
     */
    bool is_cache_outdated_;
};

} // namespace cldes

// class methods definitions
#include "cldes/src/des/DESystemCore.hpp"

// include transition proxy class
#include "cldes/TransitionProxy.hpp"

// Matrix proxy for sync operation
#include "cldes/operations/SyncSysProxy.hpp"

#endif // DESYSTEM_HPP
