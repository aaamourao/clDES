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
#include <vector>
#include "constants.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace cldes {

namespace ublas = boost::numeric::ublas;

class DESystem {
public:
    /*! \brief DESystem constructor
     *
     * Creates the DESystem object with N states defined by the argument
     * aStatesNumber and represented by its graph defined by argument aGraph.
     */
    DESystem(ublas::compressed_matrix<ScalarType> &aGraph,
             const int &aStatesNumber, const int &aInitState,
             std::vector<int> aMarkedStates, const bool &aDevCacheEnabled);

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
    ublas::compressed_matrix<ScalarType> GetGraph() const;

    /*! \brief DESystem::accessible_part() method
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a vector with all nodes
     * which are accessible from the initial state.
     */
    std::set<int> AccessiblePart();

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
    const ublas::compressed_matrix<ScalarType> *graph_;

    /*! \brief DESystem::states_number_ data member
     *
     * Hold the number of states that the automata contains. As the automata can
     * be cut, the states number is not a constant at all.
     */
    int states_number_;

    /*! \brief DESystem::dev_cache_enabled_ data member
     *
     * If dev_cache_enabled_ is true, the graph should be cached on the device
     * memory, so device_graph_ is not nullptr. It can be set at any time at run
     * time, so it is not a constant.
     */
    bool dev_cache_enabled_;

    /*! \brief DESystem::device_graph_ data member
     *
     * Transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     * TODO: Should it be a smart pointer?
     */
    viennacl::compressed_matrix<ScalarType> *device_graph_;

    /*! \brief DESystem::init_state_ data member
     *
     * Hold the initial state position.
     */
    const int init_state_;

    /*! \brief DESystem::marked_states_ data member
     *
     * Hold all marked states. Cannot be const, since the automata can be cut,
     * and some marked states may be deleted.
     */
    std::vector<int> marked_states_;

    /*! \brief DESystem::CacheGraph_ private method
     *
     * Put graph transposed data on the device memory.
     */
    void CacheGraph_();
};

}  // namespace cldes

#endif  // DESYSTEM_HPP
