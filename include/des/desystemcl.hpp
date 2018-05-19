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
 Description: DESystemCL class definition. DESystemCL is a graph, which is
 modeled as a Sparce Adjacency Matrix.
 =========================================================================
*/

#ifndef DESYSTEM_HPP
#define DESYSTEM_HPP

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <set>
#include "constants.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace cldes {

namespace ublas = boost::numeric::ublas;

/*
 * Forward declarion of DESystemCL class necessary for the forward declaration
 * of the DESystemCL's friend function op::Synchronize
 */
class DESystemCL;

namespace op {
/*
 * Forward declarion of DESystemCL's friend function Synchronize which
 * implements the parallel composition between two DES.
 */
cldes::DESystemCL Synchronize(DESystemCL &aSys0, DESystemCL &aSys1);

struct StatesTable;

StatesTable *SynchronizeStage1(DESystemCL const &aSys0,
                               DESystemCL const &aSys1);

cldes::DESystemCL SynchronizeStage2(StatesTable const *aTable,
                                    cldes::DESystemCL &aSys0,
                                    cldes::DESystemCL &aSys1);
}  // namespace op

class DESystemCL {
public:
    using GraphHostData = ublas::compressed_matrix<ScalarType>;
    using GraphDeviceData = viennacl::compressed_matrix<ScalarType>;
    using StatesSet = std::set<cldes_size_t>;
    using StatesVector = ublas::compressed_matrix<ScalarType>;
    using StatesDeviceVector = viennacl::compressed_matrix<ScalarType>;
    using StatesDenseVector = ublas::vector<ScalarType>;
    using StatesDeviceDenseVector = viennacl::vector<ScalarType>;
    using EventsSet = std::set<ScalarType>;

    /*! \brief DESystemCL constructor by copying ublas object
     *
     * Creates the DESystemCL object with N states defined by the argument
     * aStatesNumber and represented by its graph defined by argument the
     * ublas compressed matrix aGraph.
     *
     * @param aGraph Ublas matrix containing the graph data
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     * @param aMarkedStates System's marked states
     * @aDevCacheEnabled Enable or disable device cache for graph data
     */
    explicit DESystemCL(GraphHostData const &aGraph,
                        cldes_size_t const &aStatesNumber,
                        cldes_size_t const &aInitState,
                        StatesSet &aMarkedStates,
                        bool const &aDevCacheEnabled = true);

    /*! \brief DESystemCL constructor with empty matrix
     *
     * Overloads DESystemCL constructor: does not require to create a
     * ublas::compressed_matrix by the class user.
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     * @param aMarkedStates System's marked states
     * @aDevCacheEnabled Enable or disable device cache for graph data
     */
    explicit DESystemCL(cldes_size_t const &aStatesNumber,
                        cldes_size_t const &aInitState,
                        StatesSet &aMarkedStates,
                        bool const &aDevCacheEnabled = true);

    DESystemCL(DESystemCL const &aSys);

    /*! \brief DESystemCL destructor
     *
     * Delete dinamically allocated data: graph and device_graph.
     */
    virtual ~DESystemCL();

    /*! \brief Graph getter
     *
     * Returns a copy of DESystemCL's private data member graph. Considering
     * that graph is a pointer, it returns the contents of graph.
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
    DESystemCL Trim(bool const &aDevCacheEnabled = true);

    /*! \brief Returns number of states of the system
     *
     * Returns states_value_ by value.
     */
    cldes_size_t Size() const { return states_number_; }

    /*! \brief Set events_
     *
     * Set the member events_ with a set containing all events that are present
     * on the current system.
     */
    void InsertEvents(EventsSet const &aEvents);

    /*
     * TODO:
     * getters
     * enable dev cache
     * ...
     */
protected:
    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    DESystemCL();

private:
    friend DESystemCL op::Synchronize(DESystemCL &aSys0, DESystemCL &aSys1);
    friend op::StatesTable *op::SynchronizeStage1(DESystemCL const &aSys0,
                                                  DESystemCL const &aSys1);
    friend DESystemCL op::SynchronizeStage2(op::StatesTable const *aTable,
                                            DESystemCL &aSys0,
                                            DESystemCL &aSys1);

    /*! \brief Graph represented by an adjascency matrix
     *
     * A sparse matrix who represents the automata as a graph in an adjascency
     * matrix. It is implemented as a CSR scheme. The pointer is constant, but
     * its content should not be constant, as the graph should change many times
     * at runtime.
     *
     * TODO: Explain transition scheme.
     * TODO: Should it be a smart pointer?
     */
    GraphHostData *const graph_;

    /*! \brief Graph data on device memory
     *
     * Transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     * TODO: Should it be a smart pointer?
     */
    GraphDeviceData *device_graph_;

    /*! \brief Keeps if caching graph data on device is enabled
     *
     * If dev_cache_enabled_ is true, the graph should be cached on the device
     * memory, so device_graph_ is not nullptr. It can be set at any time at run
     * time, so it is not a constant.
     */
    bool dev_cache_enabled_;

    /*! \brief Keeps track if the device graph cache is outdated
     *
     * Tracks if cache, dev_graph_, needs to be updated or not.
     */
    bool is_cache_outdated_;

    /*! \brief Current system's states number
     *
     * Hold the number of states that the automata contains. As the automata can
     * be cut, the states number is not a constant at all.
     */
    cldes_size_t states_number_;

    /*! \brief Current system's initial state
     *
     * Hold the initial state position.
     */
    cldes_size_t const init_state_;

    /*! \brief Current system's marked states
     *
     * Hold all marked states. Cannot be const, since the automata can be cut,
     * and some marked states may be deleted.
     */
    StatesSet marked_states_;

    /*! \brief System's events
     *
     * A std::set containing all the events that matter for the current system.
     */
    EventsSet events_;

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
     */
    template <class StatesType>
    StatesSet *Bfs_(StatesType const &aInitialNodes,
                    std::function<void(cldes_size_t const &,
                                       cldes_size_t const &)> const &aBfsVisit);

    /*! \brief Setup BFS using dense vector and return accessed states array
     *
     * Executes a breadth first search on the graph starting from the node
     * aInitialNode. The algorithm is based on SpMV.
     *
     * @param aInitialNode Node where the searches will start
     */
    StatesSet *BfsSpMV_(
        cldes_size_t const &aInitialNode,
        std::function<void(cldes_size_t const &, cldes_size_t const &)> const
            &aBfsVisit);

    /*! \brief Calculates Bfs and returns accessed states array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpGEMM.
     *
     * @param aInitialNode Where the search will start
     */
    StatesSet *BfsCalc_(
        StatesVector &aHostX,
        std::function<void(cldes_size_t const &, cldes_size_t const &)> const
            &aBfsVisit,
        std::vector<cldes_size_t> const *const aStatesMap);

    /*! \brief Calculates Bfs using dense vector and returns accessed states
     * array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpMV.
     *
     * @param aInitialNode Where the search will start
     */
    StatesSet *BfsCalcSpMV_(
        StatesDenseVector &aHostX,
        std::function<void(cldes_size_t const &, cldes_size_t const &)> const
            &aBfsVisit);

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    StatesSet *Bfs_();
};

}  // namespace cldes
#endif  // DESYSTEM_HPP
