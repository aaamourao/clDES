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

#include "cldes/DESystemBase.hpp"
#include "cldes/src/des/DESystemCLCore.hpp"
#include "cldes/DESystem.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace cldes {

namespace op {
namespace backend {
class OclBackend;
}

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
    inline DESystemCL(DESystem const& aSys)
      : DESystem{ aSys.Size(), aSys.GetMarkedStates() }
    {
        graph_ = aSys.bit_graph_.cast<float>();
    }

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
    inline std::shared_ptr<DESystemBase> Clone() const override
    {
        std::shared_ptr<DESystemBase> this_ptr =
          std::make_shared<DESystemCL>(*this);
        return this_ptr;
    }

    /*! \brief Check if this system is a virtual proxy
     * \details DESystemCL is always a real sys
     */
    inline bool IsVirtual() const override { return false; }

    /*! \brief Graph getter
     *
     * Returns a copy of DESystemCL's private data member graph. Considering
     * that graph is a pointer, it returns the contents of graph.
     */
    GraphDeviceData GetGraph() const { return graph_ };

    /*! \brief Returns state set containing the accessible part of automa
     *
     * Executes a Breadth First Search in the graph, which represents the DES,
     * starting from its initial state. It returns a set containing all nodes
     * which are accessible from the initial state.
     */
    StatesSet AccessiblePart() const
    {
        auto accessible_states = Bfs_();
        return *accessible_states;
    }

    /*! \brief Returns state set containing the coaccessible part of automata
     *
     * Executes a Breadth First Search in the graph, until it reaches a marked
     * state.
     */
    StatesSet CoaccessiblePart()
    {
        auto coaccessible_part = Bfs_();
        return *coaccessible_part;
    }

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
    static backend::OclBackend* backend_ptr_;

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
        aBfsVisit)
    {
        /*
         * BFS on a Linear Algebra approach:
         *     Y = G^T * X
         */
        StatesVector host_x{ states_number_, 1 };

        // GPUs does not allow dynamic memory allocation. So, we have
        // to set X on host first.
        host_x(aInitialNode, 0) = 1;

        return BfsCalc_(host_x, aBfsVisit, nullptr);
    }

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
      std::vector<cldes_size_t> const* const aStatesMap)
    {
        cl_uint n_initial_nodes = aHostX.cols();

        // Copy search vector to device memory
        StatesDeviceVector x;
        viennacl::copy(aHostX, x);

        // Copy to device memory
        viennacl::copy(trans(graph_), device_graph_);

        // Executes BFS
        StatesDeviceVector y;
        auto n_accessed_states = 0;
        for (auto i = 0; i < states_number_; ++i) {
            // Using auto bellow results in compile error
            // on the following for statement
            y = viennacl::linalg::prod(dev_graph, x);

            if (n_accessed_states == y.nnz()) {
                break;
            } else {
                n_accessed_states = y.nnz();
            }

            x = y;
        }

        viennacl::copy(y, aHostX);

        // Add results to a std::set vector
        auto accessed_states = new StatesSet[n_initial_nodes];
        for (auto node = aHostX.begin1(); node != aHostX.end1(); ++node) {
            for (auto elem = node.begin(); elem != node.end(); ++elem) {
                accessed_states[elem.index2()].emplace(node.index1());
            }
        }

        return accessed_states;
    }

    /*! \brief Return a pointer to accessed states from the initial state
     *
     * Executes a breadth first search on the graph starting from init_state_.
     */
    std::shared_ptr<StatesSet> Bfs_() { return Bfs_(init_state_, nullptr); }
};

// OclBackend disable by now
backend::OclBackend* DESystemCL::backend_ptr_ = nullptr;

} // namespace cldes
#endif // DESYSTEMCL_HPP
