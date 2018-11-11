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

 File: cldes/DESystemCL.hpp
 Description: DESystemCL class definition. DESystemCL is a Discrete-Event
 System on the device memory.
 =========================================================================
*/

#ifndef DESYSTEMCL_HPP
#define DESYSTEMCL_HPP

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
template<uint8_t NEvents, typename StorageIndex = unsigned int>
class DESystemCL : public DESystemBase<NEvents, StorageIndex>
{
public:
    /*! Unfortunatelly, viennacl only accepts float type
     */
    using GraphHostData = Eigen::SparseMatrix<float, Eigen::RowMajor>;
    using StatesVector = Eigen::SparseMatrix<float, Eigen::RowMajor>;

    /*! \brief StorageIndex signed type
     * \details Eigen uses signed indexes
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Base alias
     * \details Alias to base class with implicit template params
     */
    using DESystemBase = DESystemBase<NEvents, StorageIndex>;

    /*! \brief Set of states type
     *  \details Set containg unsigned interget types which represent states.
     */
    using StatesSet = typename DESystemBase::StatesSet;

    /*! Adjacency matrix on device memory
     */
    using GraphDeviceData = viennacl::compressed_matrix<float>;

    /*! Sparse vector for BFS operations
     */
    using StatesDeviceVector = viennacl::compressed_matrix<float>;

    /*! \brief Graph const iterator
     * \details Used to iterate over the bfs result
     */
    using ColIteratorConst = Eigen::InnerIterator<StatesVector const>;

    /*! \brief DESystemCL constructor by copying ublas object
     *
     * Creates a new system based on a existent host system.
     *
     * @param aSys System on device memory
     */
    DESystemCL(DESystem<NEvents, StorageIndex>& aSys);

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

    /*! \brief Check if this system is a proxy
     * \details DESystemCL is always a real sys
     */
    bool IsVirtual() const override;

    /*! \brief Graph getter
     *
     * Returns a copy of DESystemCL's private data member graph. Considering
     * that graph is a pointer, it returns the contents of graph.
     */
    GraphDeviceData getGraph() const;

    /*! \brief Returns state set containing the accessible part of automa
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSet accessiblePart();

    /*! \brief Returns state set containing the coaccessible part of automata
     *
     * Executes a Breadth First Search in the graph, until it reaches a marked
     * state.
     */
    StatesSet coaccessiblePart();

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool Containstrans(StorageIndex const& aQ,
                       ScalarType const& aEvent) const override;

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StorageIndexSigned trans(StorageIndex const& aQ,
                             ScalarType const& aEvent) const override;

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool Containsinvtrans(StorageIndex const& aQ,
                          ScalarType const& aEvent) const override;

    /*! \brief Returns DES inverse transition: q = f^-1(q_to, e)
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    StatesArray<StorageIndex> invtrans(StorageIndex const& aQfrom,
                                       ScalarType const& aEvent) const override;

    /*! \brief Returns EventsSet relative to state q
     *
     * @param aQ A state on the sys
     */
    EventsSet<NEvents> GetStateEvents(StorageIndex const& aQ) const override;

    /*! \brief Returns EventsSet relative to state inv q
     *
     * @param aQ A state on the sys
     */
    EventsSet<NEvents> GetInvStateEvents(StorageIndex const& aQ) const override;

    /*! \brief Invert graph
     *
     */
    void AllocateInvertedGraph() const override;

    /*! \brief Free inverted graph
     *
     */
    void ClearInvertedGraph() const override;

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
     * transposed graph_ data, but on device memory (usually a GPU). It is a
     * dev_cache_enabled_ is false. It cannot be const, since it may change as
     * dev_cache_enabled_ changes.
     *
     */
    GraphDeviceData device_graph_;
    GraphHostData graph_;

    /*! \brief Pointer to beckend singleton object
     */
    static std::unique_ptr<backend::OclBackend> backend_ptr_;

    /*! \brief Setup BFS and return accessed states array
     *
     * Executes a breadth first search on the graph starting from N nodes
     * in aInitialNodes. The algorithm is based on SpGEMM.
     *
     * @param aInitialNodes Set of nodes where the searches will start
     */
    template<class StatesType>
    std::shared_ptr<StatesSet> bfs_(
      StatesType const& aInitialNodes,
      std::function<void(StorageIndex const&, StorageIndex const&)> const&
        aBfsVisit);

    /*! \brief Calculates Bfs and returns accessed states array
     *
     * Executes a breadth first search on the graph starting from one single
     * node. The algorithm is based on SpGEMM.
     *
     * @param aInitialNode Where the search will start
     */
    std::shared_ptr<StatesSet> bfsCalc_(
      StatesVector& aHostX,
      std::function<void(StorageIndex const&, StorageIndex const&)> const&
        aBfsVisit,
      std::vector<StorageIndex> const* const aStatesMap);

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    std::shared_ptr<StatesSet> bfs_();
};

} // namespace cldes

#include "cldes/src/des/DESystemCLCore.hpp"

#endif // DESYSTEMCL_HPP
