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

    /*! \brief DESystem constructor
     *
     * Creates the DESystem object with N states defined by the argument
     * aStatesNumber and represented by its graph defined by argument aGraph.
     */
    DESystem(GraphHostData &aGraph, cldes_size_t const &aStatesNumber,
             cldes_size_t const &aInitState, StatesSet &aMarkedStates,
             bool const &aDevCacheEnabled = true);

    /*! \brief DESystem constructor
     *
     * Overloads DESystem constructor: does not require to create a
     * ublas::compressed_matrix by the class user.
     */
    DESystem(cldes_size_t const &aStatesNumber, cldes_size_t const &aInitState,
             StatesSet &aMarkedStates, bool const &aDevCacheEnabled = true);

    /*! \brief DESystem destructor
     *
     * Delete dinamically allocated data: graph and device_graph.
     */
    virtual ~DESystem();

    /*! \brief DESystem::getgraph() method
     *
     * Returns a copy of DESystem's private data member graph. Considering that
     * graph is a pointer, it returns the contents of graph.
     */
    GraphHostData GetGraph() const;

    /*! \brief DESystem::accessible_part() method
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSet AccessiblePart();

    /*! \brief DESystem::operator()
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

private:
    /*! \brief DESystem::graph_ data member
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

    /*! \brief DESystem::device_graph_ data member
     *
     * Transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     * TODO: Should it be a smart pointer?
     */
    GraphDeviceData *device_graph_;

    /*! \brief Device Cache Enabled private data member
     *
     * If dev_cache_enabled_ is true, the graph should be cached on the device
     * memory, so device_graph_ is not nullptr. It can be set at any time at run
     * time, so it is not a constant.
     */
    bool dev_cache_enabled_;

    /*! \brief Is Cache Outdated private data member
     *
     * Tracks if cache, dev_graph_, needs to be updated or not.
     */
    bool is_cache_outdated_;

    /*! \brief States Number private data member
     *
     * Hold the number of states that the automata contains. As the automata can
     * be cut, the states number is not a constant at all.
     */
    cldes_size_t states_number_;

    /*! \brief Initial State private data member
     *
     * Hold the initial state position.
     */
    cldes_size_t const init_state_;

    /*! \brief Marked States private data member
     *
     * Hold all marked states. Cannot be const, since the automata can be cut,
     * and some marked states may be deleted.
     */
    StatesSet marked_states_;

    /*! \brief Cache Graph private method
     *
     * Put graph transposed data on the device memory.
     */
    void CacheGraph_();

    /*! \brief Upgrade Graph Cache private method
     *
     * Refresh the graph data on device memory.
     */
    void UpdateGraphCache_();
};

}  // namespace cldes

#endif  // DESYSTEM_HPP
