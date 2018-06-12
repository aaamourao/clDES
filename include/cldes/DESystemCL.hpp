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

 File: cldes/DESystemCL.hpp
 Description: DESystemCL class definition. DESystemCL is a Discrete-Event
 System on the device memory.
 =========================================================================
*/

#ifndef DESYSTEMCL_HPP
#define DESYSTEMCL_HPP

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include "cldes/Constants.hpp"
#include "cldes/DESystemCLCore.hpp"
#include "viennacl/compressed_matrix.hpp"
#include <Eigen/Sparse>
#include <set>

namespace cldes {

namespace op {
struct StatesTable;

/*
 * Forward declarion of DESystemCL's friend function Synchronize which
 * implements the parallel composition between two DES.
 */
cldes::DESystemCL
Synchronize(DESystemCL& aSys0, DESystemCL& aSys1);

StatesTable*
SynchronizeStage1(DESystemCL const& aSys0, DESystemCL const& aSys1);

cldes::DESystemCL
SynchronizeStage2(StatesTable const* aTable,
                  cldes::DESystemCL& aSys0,
                  cldes::DESystemCL& aSys1);
} // namespace op

namespace backend {
class OclBackend;
}

/*! \brief Discrete-Events System on device memory
 *
 * Implement a DES on the device memory and their respective operations for GPUs.
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystemCL
{
public:
    using GraphDeviceData = viennacl::compressed_matrix<ScalarType>;
    using StatesDeviceVector = viennacl::compressed_matrix<ScalarType>;
    using StatesSet = std::set<cldes_size_t>;

    /*! \brief DESystemCL constructor by copying ublas object
     *
     * Creates a new system based on a existent host system.
     *
     * @param aSys System on device memory
     */
    explicit DESystemCL(DESystem const& aSys);

    /*! \brief DESystemCL destructor
     */
    // virtual ~DESystemCL();

    /*! \brief Graph getter
     *
     * Returns a copy of DESystemCL's private data member graph. Considering
     * that graph is a pointer, it returns the contents of graph.
     */
    GraphDeviceData GetGraph() const;

    /*! \brief Returns state set containing the accessible part of automa
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSetDevice AccessiblePart();

    /*! \brief Returns state set containing the coaccessible part of automata
     *
     * Executes a Breadth First Search in the graph, until it reaches a marked
     * state.
     */
    StatesSetDevice CoaccessiblePart();

    /*! \brief Returns States Set which is the Trim part of the system
     *
     * Gets the intersection between the accessible part and the coaccessible
     * part.
     */
    StatesSetDevice TrimStates();

    /*! \brief Returns DES which is the Trim part of this
     *
     * Cut the non-accessible part of current system and then cut the
     * non-coaccessible part of the last result. The final resultant system
     * is called a trim system.
     */
    DESystemCL Trim();

    /*! \brief Returns number of states of the system
     *
     * Returns states_value_ by value.
     */
    StorageIndex Size() const { return states_number_; }

    /*! \brief Set events_
     *
     * Set the member events_ with a set containing all events that are present
     * on the current system.
     */
    void InsertEvents(EventsSet const& aEvents);

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
    friend DESystemCL op::Synchronize(DESystemCL& aSys0, DESystemCL& aSys1);
    friend op::StatesTable* op::SynchronizeStage1(DESystemCL const& aSys0,
                                                  DESystemCL const& aSys1);
    friend DESystemCL op::SynchronizeStage2(op::StatesTable const* aTable,
                                            DESystemCL& aSys0,
                                            DESystemCL& aSys1);

    /*! \brief Graph data on device memory
     *
     * Transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     */
    GraphDeviceData device_graph_;

    /*! \brief Current system's states number
     *
     * Hold the number of states that the automata contains. As the automata can
     * be cut, the states number is not a constant at all.
     */
    StorageIndex states_number_;

    /*! \brief Current system's initial state
     *
     * Hold the initial state position.
     */
    StorageIndex init_state_;

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

    static backend::OclBackend* backend_ptr_;

    /*! \brief Setup BFS and return accessed states array
     *
     * Executes a breadth first search on the graph starting from N nodes
     * in aInitialNodes. The algorithm is based on SpGEMM.
     *
     * @param aInitialNodes Set of nodes where the searches will start
     */
    template<class StatesType>
    StatesSet* Bfs_(StatesType const& aInitialNodes,
                    std::function<void(cldes_size_t const&,
                                       cldes_size_t const&)> const& aBfsVisit);

    /*! \brief Setup BFS using dense vector and return accessed states array
     *
     * Executes a breadth first search on the graph starting from the node
     * aInitialNode. The algorithm is based on SpMV.
     *
     * @param aInitialNode Node where the searches will start
     */
    StatesSet* BfsSpMV_(
      cldes_size_t const& aInitialNode,
      std::function<void(cldes_size_t const&, cldes_size_t const&)> const&
        aBfsVisit);

    /*! \brief Calculates Bfs and returns accessed states array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpGEMM.
     *
     * @param aInitialNode Where the search will start
     */
    StatesSet* BfsCalc_(
      StatesVector& aHostX,
      std::function<void(cldes_size_t const&, cldes_size_t const&)> const&
        aBfsVisit,
      std::vector<cldes_size_t> const* const aStatesMap);

    /*! \brief Calculates Bfs using dense vector and returns accessed states
     * array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpMV.
     *
     * @param aInitialNode Where the search will start
     */
    StatesSet* BfsCalcSpMV_(
      StatesDenseVector& aHostX,
      std::function<void(cldes_size_t const&, cldes_size_t const&)> const&
        aBfsVisit);

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    StatesSet* Bfs_();
};

} // namespace cldes
#endif // DESYSTEMCL_HPP
