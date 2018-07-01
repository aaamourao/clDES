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

#include "cldes/DESystem.hpp"
#include "cldes/DESystemBase.hpp"
#include "cldes/src/des/DESystemCLFwd.hpp"

namespace cldes {

/*! \class DESystemCL
 * \brief Discrete-Events System on device memory
 *
 * Implement a DES on the device memory and their respective operations for
 * GPUs.
 *
 * \tparam NEvents Number of events
 * \tparam StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystemCL : public DESystemBase<NEvents, StorageIndex>
{
public:
    /*! Unfortunatelly, viennacl only accepts float type
     */
    using GraphHostData = Eigen::SparseMatrix<float, Eigen::RowMajor>;
    using StatesVector = Eigen::SparseMatrix<float, Eigen::RowMajor>;

    /*! \brief Base alias
     * \details Alias to base class with implicit template params
     */
    using DESystemBase = DESystemBase<NEvents, StorageIndex>;
    using DESystem = DESystem<NEvents, StorageIndex>;

    /*! Adjacency matrix on device memory
     */
    using GraphDeviceData = viennacl::compressed_matrix<ScalarType>;

    /*! Sparse vector for BFS operations
     */
    using StatesDeviceVector = viennacl::compressed_matrix<ScalarType>;

    /*! \brief DESystemCL constructor by copying ublas object
     *
     * Creates a new system based on a existent host system.
     *
     * @param aSys System on device memory
     */
    DESystemCL(DESystem const& aSys);

    /*! \brief DESystemCL destructor
     */
    ~DESystemCL() = default;

    DESystemCL(DESystemCL&&) = default;

    DESystemCL(DESystemCL const&) = default;

    DESystemCL<NEvents, StorageIndex>& operator=(DESystemCL&&) = default;

    DESystemCL<NEvents, StorageIndex>& operator=(DESystemCL const&) = default;

    /*! \brief Clone method for polymorphic copy
     *  \return Shared pointer to this object
     */
    std::shared_ptr<DESystemBase> Clone() const override;

    /*! \brief Check if this system is a virtual proxy
     * \details DESystemCL is always a real sys
     */
    bool IsVirtual() const override;

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
    StatesSet AccessiblePart() const;

    /*! \brief Returns state set containing the coaccessible part of automata
     *
     * Executes a Breadth First Search in the graph, until it reaches a marked
     * state.
     */
    StatesSet CoaccessiblePart() const;

protected:
    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    explicit DESystemCL() = default;

private:
    friend class DESystem<NEvents, StorageIndex>;

    /*! \brief Graph data on device memory
     *
     * Transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     */
    GraphDeviceData device_graph_;
    GraphHostData graph_;

    /*! \brief Pointer to beckend singleton object
     */
    static std::shared_ptr<backend::OclBackend> backend_ptr_;

    /*! \brief Setup BFS and return accessed states array
     *
     * Executes a breadth first search on the graph starting from N nodes
     * in aInitialNodes. The algorithm is based on SpGEMM.
     *
     * @param aInitialNodes Set of nodes where the searches will start
     */
    template<class StatesType>
    std::shared_ptr<StatesSet> Bfs_(
      StatesType const& aInitialNodes,
      std::function<void(cldes_size_t const&, cldes_size_t const&)> const&
        aBfsVisit);

    /*! \brief Calculates Bfs and returns accessed states array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpGEMM.
     *
     * @param aInitialNode Where the search will start
     */
    std::shared_ptr<StatesSet> BfsCalc_(
      StatesVector& aHostX,
      std::function<void(cldes_size_t const&, cldes_size_t const&)> const&
        aBfsVisit,
      std::vector<cldes_size_t> const* const aStatesMap);

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    std::shared_ptr<StatesSet> Bfs_();
};

} // namespace cldes

#include "cldes/src/des/DESystemCLCore.hpp"

#endif // DESYSTEMCL_HPP
