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

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <set>
#include "constants.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace cldes {

namespace ublas = boost::numeric::ublas;

class DESystem {
public:
    using GraphHostData = ublas::compressed_matrix<ScalarType>;
    using GraphDeviceData = viennacl::compressed_matrix<ScalarType>;
    using StatesSet = std::set<cldes_size_t>;
    using StatesVector = ublas::compressed_matrix<ScalarType>;
    using StatesDeviceVector = viennacl::compressed_matrix<ScalarType>;

    /*! \brief DESystem constructor by copying ublas object
     *
     * Creates the DESystem object with N states defined by the argument
     * aStatesNumber and represented by its graph defined by argument the
     * ublas compressed matrix aGraph.
     */
    explicit DESystem(GraphHostData const &aGraph,
                      cldes_size_t const &aStatesNumber,
                      cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                      bool const &aDevCacheEnabled = true);

    /*! \brief DESystem constructor with empty matrix
     *
     * Overloads DESystem constructor: does not require to create a
     * ublas::compressed_matrix by the class user.
     */
    explicit DESystem(cldes_size_t const &aStatesNumber,
                      cldes_size_t const &aInitState, StatesSet &aMarkedStates,
                      bool const &aDevCacheEnabled = true);

    /*! \brief DESystem destructor
     *
     * Delete dinamically allocated data: graph and device_graph.
     */
    virtual ~DESystem();

    /*! \brief Graph getter
     *
     * Returns a copy of DESystem's private data member graph. Considering that
     * graph is a pointer, it returns the contents of graph.
     */
    GraphHostData GetGraph() const;

    /*! \brief Automata accessible part operation
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSet AccessiblePart();

    /*! \brief Operator "()" for assigning values to elements
     *
     * Override operator () for changing transinstions with a single assignment:
     * e.g. discrete_system_foo(2,1) = 3.0f;
     */
    GraphHostData::reference operator()(cldes_size_t const &lin,
                                        cldes_size_t const &col);

    /*
     * TODO:
     * getters
     * enable dev cache
     * CoaccessiblePart
     * Trim
     * ...
     */
protected:
    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    DESystem();

private:
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

    /*! \brief Breadth First Search
     *
     * Executes a breadth first search on the graph.
     */
    StatesSet *Bfs_(cldes_size_t const &aInitialNode);

    /*! \brief Breadth First Search
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    StatesSet *Bfs_() { Bfs_(init_state_); };
};

}  // namespace cldes

#endif  // DESYSTEM_HPP
