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

 File: desystem.hpp
 Description: DESystem class definition. DESystem is a graph, which is
 modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#ifndef DESYSTEM_HPP
#define DESYSTEM_HPP

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include "cldes/Constants.hpp"
#include <Eigen/Sparse>
#include <set>
#include <sparsepp/spp.h>
#include <stack>
#include <tuple>
#include <vector>

// Forward declarations
namespace cldes {
/*
 * Forward declarion of DESystem class necessary for the forward declaration of
 * the DESystem's friend function op::Synchronize
 */
template<size_t NEvents = 32u, typename StorageIndex = unsigned>
class DESystem;

/*
 * Forward declarion of DESystem's friends class TransitionProxy. A transition
 * is an element of the adjascency matrix which implements the des graph.
 */
template<size_t NEvents, typename StorageIndex>
class TransitionProxy;

// Forward declartions of friends functions which implement des operations
namespace op {
template<size_t NEvents, typename StorageIndex>
using GraphType =
  Eigen::SparseMatrix<std::bitset<NEvents>,
                      Eigen::RowMajor,
                      typename std::make_signed<StorageIndex>::type>;

/*
 * Forward declarion of DESystem's friend function Synchronize which
 * implements the parallel composition between two DES.
 */
template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
Synchronize(DESystem<NEvents, StorageIndex> const& aSys0,
            DESystem<NEvents, StorageIndex> const& aSys1);

// Table containing states of a system on the device mem
struct StatesTable;

/*
 * Tuple which implements the virtual states of sync a system on the device
 * memory
 */
struct StatesTuple;

// Tuple representing a virtual state
template<typename StorageIndex>
using StatesTupleHost = std::pair<StorageIndex, StorageIndex>;

// Table containing virtual states of a virtual system as integers/longs
template<typename StorageIndex>
using StatesTableHost = spp::sparse_hash_set<StorageIndex>;

// Stack of states for implementing a DFS
template<typename StorageIndex>
using StatesStack = std::stack<StorageIndex>;

// Table of events
using EventsTableHost = spp::sparse_hash_set<unsigned>;

// Implements the first stage of the parallel composition: gen tuples
template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SynchronizeStage1(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

// Implements the second stage of the parallel composition: calc transitions
template<size_t NEvents, typename StorageIndex>
void
SynchronizeStage2(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1);

// Calculate f(q, event) of a virtual system
template<size_t NEvents, typename StorageIndex>
StorageIndex
TransitionVirtual(cldes::DESystem<NEvents, StorageIndex> const& aSys0,
                  cldes::DESystem<NEvents, StorageIndex> const& aSys1,
                  StorageIndex const& aQ,
                  cldes::ScalarType const& aEvent);

// Remove bad states recursively
template<size_t NEvents, typename StorageIndex>
void
RemoveBadStates(cldes::DESystem<NEvents, StorageIndex>& aVirtualSys,
                cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aE,
                GraphType<NEvents, StorageIndex> const& aInvGraphP,
                GraphType<NEvents, StorageIndex> const& aInvGraphE,
                StatesTableHost<StorageIndex>& aC,
                StorageIndex const& aQ,
                std::bitset<NEvents> const& aNonContrBit,
                StatesTableHost<StorageIndex>& aRmTable);

// Calculate the monolithic supervisor of a fiven plant and spec
template<size_t NEvents, typename StorageIndex>
cldes::DESystem<NEvents, StorageIndex>
SupervisorSynth(cldes::DESystem<NEvents, StorageIndex> const& aP,
                cldes::DESystem<NEvents, StorageIndex> const& aS,
                EventsTableHost const& aNonContr);
} // namespace op

/*! \brief Discrete-Events System on host memory
 *
 * Implement a DES on the host memory and their respective operations for CPUs.
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<size_t NEvents, typename StorageIndex>
class DESystem
{
public:
    /*! \brief Bit array representing an events set
     *
     * Each bit represent a different event.
     * 0 -> does not contain event
     * 1 -> contains event
     */
    using EventsSet = std::bitset<NEvents>;

    /*! \brief Adjacency matrix of bitarrays implementing a graph
     *
     * The graph represents the DES automata:
     * Non zero element: contains at least a transition
     * Each non 0 bit of each element: event that lead to the next stage
     *
     * row index: from state
     * col index: to state
     */
    using GraphHostData =
      Eigen::SparseMatrix<EventsSet,
                          Eigen::RowMajor,
                          typename std::make_signed<StorageIndex>::type>;

    /*! \brief Adjacency matrix of bit implementing a graph
     *
     * Sparse matrix implementing a graph with an adjacency matrix.
     * The graph represents a simplified DES automata:
     * Non zero element (true): contains at least a transition
     *
     * row index: from state
     * col index: to state
     */
    using BitGraphHostData =
      Eigen::SparseMatrix<bool,
                          Eigen::RowMajor,
                          typename std::make_signed<StorageIndex>::type>;

    /*! \brief Adjacency  matrix of bit implementing searching nodes
     *
     * Structure used for traversing the graph using a linear algebra approach
     */
    using StatesVector =
      Eigen::SparseMatrix<bool,
                          Eigen::ColMajor,
                          typename std::make_signed<StorageIndex>::type>;

    /*! \brief Set of states type
     */
    using StatesSet = std::set<StorageIndex>;

    /*! \brief Table of transitions on a STL container
     */
    using StatesEventsTable = std::vector<EventsSet>;

    /*! \brief Set of Events implemented as a Hash Table for searching
     * efficiently.
     */
    using EventsTable = spp::sparse_hash_set<unsigned>;

    /*! \brief Vector of states type
     */
    using StatesTable = std::vector<StorageIndex>;

    /*! \brief Arguments of a transition function vector
     *
     * f(s, e) = s_out -> (s, e) are the arguments
     */
    using ArgTransition = std::vector<std::pair<StorageIndex, EventsSet>> ;

    /*! \brief Vector of inverted transitions
     *
     * f(s, e) = s_out -> (s_out, (s, e)) is the inverted transition.
     */
    using InvTransition = std::vector<std::pair<StorageIndex, ArgTransition>>;

    /*! \brief DESystem constructor with empty matrix
     *
     * Overloads DESystem constructor: does not require to create a
     * eigen sparse matrix by the class user.
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     * @param aMarkedStates System's marked states
     * @aDevCacheEnabled Enable or disable device cache for graph data
     */
    explicit DESystem(StorageIndex const& aStatesNumber,
                      StorageIndex const& aInitState,
                      StatesSet& aMarkedStates,
                      bool const& aDevCacheEnabled = true);

    /*! \brief Copy constructor
     *
     * TODO: It will be necessary when the members used only for virtual systems
     * become pointers.
     */
    // DESystem(DESystem &aSys);

    /*! \brief DESystem destructor
     *
     * TODO: It will be necessary when the members used only for virtual systems
     * become pointers.
     */
    // virtual ~DESystem();

    /*! \brief Graph getter
     *
     * Returns a copy of DESystem's private data member graph. Considering that
     * graph is a pointer, it returns the contents of graph.
     */
    GraphHostData GetGraph() const;

    /*! \brief Returns state set containing the accessible part of automa
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSet AccessiblePart();

    /*! \brief Returns state set containing the coaccessible part of automata
     *
     * Executes a Breadth First Search in the graph, until it reaches a marked
     * state.
     */
    StatesSet CoaccessiblePart();

    /*! \brief Returns States Set which is the Trim part of the system
     *
     * Gets the intersection between the accessible part and the coaccessible
     * part.
     */
    StatesSet TrimStates();

    /*! \brief Returns DES which is the Trim part of this
     *
     * Cut the non-accessible part of current system and then cut the
     * non-coaccessible part of the last result. The final resultant system
     * is called a trim system.
     *
     * @param aDevCacheEnabled Enables cache device graph on returned DES
     */
    void Trim();

    /*! \brief Returns value of the specified transition
     *
     * Override operator () for reading transinstions values:
     * e.g. discrete_system_foo(2,1);
     *
     * @param aLin Element's line
     * @param aCol Element's column
     */
    EventsSet const operator()(StorageIndex const& aLin,
                               StorageIndex const& aCol) const;

    /*! \brief Returns value of the specified transition
     *
     * Override operator () for changing transinstions with a single assignment:
     * e.g. discrete_system_foo(2,1) = 3.0f;
     *
     * @param aLin Element's line
     * @param aCol Element's column
     */
    TransitionProxy<NEvents, StorageIndex> operator()(StorageIndex const& aLin,
                                                      StorageIndex const& aCol);

    /*! \brief Returns number of states of the system
     *
     * Returns states_value_ by value.
     */
    StorageIndex Size() const { return states_number_; }

    /*! \brief Insert events
     *
     * Set the member events_ with a set containing all events that are present
     * on the current system. Take care using this method. It was designed for
     * testing and debugging.
     *
     * @params aEvents Set containing all the new events of the current system
     */
    void InsertEvents(EventsSet const& aEvents);

    /*
     * TODO:
     * getters
     * setters: e.g. remove transition, remove state
     * enable dev cache
     * ...
     */
protected:
    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    DESystem(){};

private:
    /*! \brief DESystem operations
     *
     * Functions which implement DES operations between two systems and need
     * to access private members of DESystem. The design constrains that lead to
     * declare all these friends functions are:
     *
     * * Functions which have more than one system as input and returns a
     * different system should not be a member of DESystem class.
     * * More efficiency when accessing vars than using getters and setters.
     * * Define a different namespace for DES operations.
     */
    friend class TransitionProxy<NEvents, StorageIndex>;
    friend DESystem cldes::op::Synchronize<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex> const& aSys0,
      DESystem<NEvents, StorageIndex> const& aSys1);
    friend DESystem cldes::op::SynchronizeStage1<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex> const& aSys0,
      DESystem<NEvents, StorageIndex> const& aSys1);
    friend void cldes::op::SynchronizeStage2<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex>& aVirtualSys,
      DESystem<NEvents, StorageIndex> const& aSys0,
      DESystem<NEvents, StorageIndex> const& aSys1);
    friend StorageIndex cldes::op::TransitionVirtual<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex> const& aSys0,
      DESystem<NEvents, StorageIndex> const& aSys1,
      StorageIndex const& aQ,
      cldes::ScalarType const& aEvent);
    friend void cldes::op::RemoveBadStates<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex>& aVirtualSys,
      DESystem<NEvents, StorageIndex> const& aP,
      DESystem<NEvents, StorageIndex> const& aE,
      op::GraphType<NEvents, StorageIndex> const& aInvGraphP,
      op::GraphType<NEvents, StorageIndex> const& aInvGraphE,
      op::StatesTableHost<StorageIndex>& aC,
      StorageIndex const& aQ,
      std::bitset<NEvents> const& aNonContr,
      op::StatesTableHost<StorageIndex>& aRmTable);
    friend DESystem cldes::op::SupervisorSynth<NEvents, StorageIndex>(
      DESystem<NEvents, StorageIndex> const& aP,
      DESystem<NEvents, StorageIndex> const& aE,
      op::EventsTableHost const& aNonContr);

    /*! \brief Graph represented by an adjascency matrix
     *
     * A sparse matrix which represents the automata as a graph in an
     * adjascency matrix. It is implemented as a CSR scheme.
     *
     * Non zero element: transition from <row index> to <col index> when events
     * represented by the set bit indexes occurs.
     *
     * e.g. M(2, 3) = 101; Transition from state 2 to state 3 with the condition
     * event 0 OR event 2.
     */
    GraphHostData graph_;

    /*! \brief Graph represented by an adjascency matrix
     *
     * A sparse bit matrix which represents when a state has any transition to
     * other states plus the identity matrix. It is used to calculate the
     * accessible part, coaccessible part and trim operations efficiently. The
     * matrix is also transposed for making accessible part faster: when
     * calculating trim states, it is always necessary to calculate the
     * accessible part first. It implies that the accessible part usually is
     * calculated with a larger matrix.
     */
    BitGraphHostData bit_graph_;

    /*! \brief Keeps if caching graph data on device is enabled
     *
     * If dev_cache_enabled_ is true, the graph should be cached on the
     * device memory, so device_graph_ is not nullptr. It can be set at any
     * time at run time, so it is not a constant.
     */
    bool dev_cache_enabled_;

    /*! \brief Keeps track if the device graph cache is outdated
     *
     * Tracks if cache, dev_graph_, needs to be updated or not.
     */
    bool is_cache_outdated_;

    /*! \brief Current system's states number
     *
     * Hold the number of states that the automata contains. As the automata
     * can be cut, the states number is not a constant at all.
     */
    StorageIndex states_number_;

    /*! \brief Current system's initial state
     *
     * Hold the initial state position.
     */
    StorageIndex init_state_;

    /*! \brief Current system's marked states
     *
     * Hold all marked states. Cannot be const, since the automata can be
     * cut, and some marked states may be deleted.
     */
    StatesSet marked_states_;

    /*! \brief System's events hash table
     *
     * A hash table containing all the events that matter for the current
     * system.
     */
    EventsSet events_;

    /*! \brief Vector containing a events hash table per state
     */
    StatesEventsTable states_events_;

    /*! \brief Virtual states contained in the current system
     *
     * Valid onnly when this system is virtual.
     * TODO: Change it to a pointer and allocate only when the system is virtual
     * and deallocate when it become a concrete system.
     */
    StatesTable virtual_states_;

    /*! \brief Events contained only in the left operator of a synchronizing op.
     *
     * Valid onnly when this system is virtual.
     * TODO: Change it to a pointer and allocate only when the system is virtual
     * and deallocate when it become a concrete system.
     */
    EventsSet only_in_0_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     *
     * Valid onnly when this system is virtual.
     * TODO: Change it to a pointer and allocate only when the system is virtual
     * and deallocate when it become a concrete system.
     */
    EventsSet only_in_1_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     *
     * Valid onnly when this system is virtual.
     * TODO: Change it to a pointer and allocate only when the system is virtual
     * and deallocate when it become a concrete system.
     */
    InvTransition transtriplet_;

    /*! \brief Vector containing a events hash table per state
     *
     * It represents the transitions of the inverted graph for the supervisor
     * synthesis.
     */
    StatesEventsTable inv_states_events_;

    /*! \brief Method for caching the graph
     *
     * Put graph transposed data on the device memory.
     */
    void CacheGraph_();

    /*! \brief Method for updating the graph
     *
     * Refresh the graph data on device memory.
     */
    void UpdateGraphCache_();

    /*! \brief Setup BFS and return accessed states array
     *
     * Executes a breadth first search on the graph starting from N nodes
     * in aInitialNodes. The algorithm is based on SpGEMM.
     *
     * @param aInitialNodes Set of nodes where the searches will start
     * @param aBfsVisit Function to execute on every accessed state: it can be
     * null
     */
    template<class StatesType>
    StatesSet* Bfs_(StatesType const& aInitialNodes,
                    std::function<void(StorageIndex const&,
                                       StorageIndex const&)> const& aBfsVisit);

    /*! \brief Overload Bfs_ for the special case of a single initial node
     */
    StatesSet* Bfs_(StorageIndex const& aInitialNode,
                    std::function<void(StorageIndex const&,
                                       StorageIndex const&)> const& aBfsVisit);

    /*! \brief Calculates Bfs and returns accessed states array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpGEMM.
     *
     * @param aInitialNode Where the search will start
     * @param aBfsVisit Function to execute on every accessed state: it can be
     * null
     * @param aStatesMap Map each vector used to implement in a bfs when
     * multiple ones are execute togheter: it can be null
     */
    StatesSet* BfsCalc_(
      StatesVector& aHostX,
      std::function<void(StorageIndex const&, StorageIndex const&)> const&
        aBfsVisit,
      std::vector<StorageIndex> const* const aStatesMap);

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from
     * init_state_.
     */
    StatesSet* Bfs_();
};

} // namespace cldes

// Include DESystem implementation
#include "cldes/src/des/DESystemCore.hpp"

#endif // DESYSTEM_HPP
