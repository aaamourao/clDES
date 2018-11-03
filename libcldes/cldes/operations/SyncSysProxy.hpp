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

 File: cldes/operations/SyncSysProxy.hpp
 Description: Virtual Proxy for multiple parallel compositions
 =========================================================================
*/
/*!
 * \file cldes/operations/SyncSysProxy.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-06-16
 *
 * Virtual Proxy for multiple parallel compositions operations.
 */

#ifndef SYNC_SYS_PROXY_HPP
#define SYNC_SYS_PROXY_HPP

#include "cldes/Constants.hpp"
#include "cldes/DESystemBase.hpp"
#include "cldes/EventsSet.hpp"
#include "cldes/src/operations/SyncSysProxyFwd.hpp"

namespace cldes {
namespace op {

/*! \class SyncSysProxy
 * \brief Parallel composition return type
 * \details Represents a parallel composition. Allows lazy evaluation, but can
 * be converted to a concrete system any time with a DESystem() cast. It is
 * implemented as a binary tree of implicit parallel compositions, where the
 * leafs are concrete systems.
 *
 * \tparam NEvents Number of events
 * \tparam StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class SyncSysProxy
  : public DESystemBase<NEvents,
                        StorageIndex,
                        SyncSysProxy<NEvents, StorageIndex>>
{
public:
    /*! \brief Signed template parameter type for eigen indexes
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Base alias
     * \details Alias to implicit speciallization of base class.
     */
    using DESystemBase =
      DESystemBase<NEvents, StorageIndex, DESystem<NEvents, StorageIndex>>;

    /*! \brief Vector of states type
     * \details Vector containing states represented by unsigned indexes.
     */
    using StatesTable = typename DESystemBase::StatesTable;

    /*! \brief Vector of inverted transitions
     * \details Vector which stores transitions:
     * f(s, e) = s_out -> (s_out, (s, e)) is the inverted transition.
     */
    using TrVector =
      std::vector<std::pair<StorageIndex, InvArgTrans<StorageIndex>*>>;

    /*! \brief Alias to the the related DESystem template
     * \details alias Alias to implicit speciallization to concrete system.
     */
    using DESystem = DESystem<NEvents, StorageIndex>;

    /*! \brief SyncSysProxy unique constructor
     * Create a binary tree that represents multiple operations.
     *
     * @param aSys0 Left operand DESystem reference
     * @param aSys1 Right operand DESystem reference
     */
    SyncSysProxy(DESystemBase const& aSys0, DESystemBase const& aSys1);

    /*! \brief DESystem destructor
     * \details Override base destructor.
     */
    ~SyncSysProxy() = default;

    /*! \brief Move constructor
     * \details Enable move semantics.
     */
    SyncSysProxy(SyncSysProxy&&) = default;

    /*! \brief Copy constructor
     * \details Enable copy by reference.
     */
    SyncSysProxy(SyncSysProxy const&) = default;

    /*! \brief Operator =
     * Use move semantics when assigning to rvalues.
     */
    SyncSysProxy<NEvents, StorageIndex>& operator=(SyncSysProxy&&) = default;

    /*! \brief Operator = to const type
     * \details Enables copy by reference.
     */
    SyncSysProxy<NEvents, StorageIndex>& operator=(SyncSysProxy const&) =
      default;

    /*! \brief Overload conversion to DESystem
     * \details Convert the current virtual proxy to a concrete system.
     */
    operator DESystem();

    /*! \brief Convert to DESystem from const proxy
     * \details Convert it to a concrete system when the current system is
     * const.
     * \warning Expensive, because it implies a copy. Avoid it.
     */
    operator DESystem() const;

    /*! \brief Is it real?
     * \details Nooo!!!
     *
     * \return False, always.
     */
    bool constexpr static isVirtual_impl() { return true; }

    /*! \brief Clone method for polymorphic copy
     * \details Method for cloning on a polymorphic way.
     *
     * \return shared pointer to base.
     */
    std::shared_ptr<DESystemBase> clone_impl() const;

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool containsTrans_impl(StorageIndex const& aQ,
                            ScalarType const& aEvent) const;

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StorageIndexSigned trans_impl(StorageIndex const& aQ,
                                  ScalarType const& aEvent) const;

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool containsInvTrans_impl(StorageIndex const& aQ,
                               ScalarType const& aEvent) const;

    /*! \brief Returns DES inverse transition
     * \details  q = f^-1(q_to, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StatesArray<StorageIndex> invTrans_impl(StorageIndex const& aQ,
                                            ScalarType const& aEvent) const;

    /*! \brief Overload operator -> to access concrete sys
     * \warning If it was not converted yet, returns nullptr.
     *
     * \return Shared pointer to concrete system.
     */
    std::shared_ptr<DESystem> operator->() { return *sys_ptr_; }

    /*! \brief Overload operator dereference to access concrete sys
     * \warning If it was not converted yet, returns nullptr.
     *
     * \return Shared pointer to concrete system.
     */
    std::shared_ptr<DESystem> operator*() { return *sys_ptr_; }

    /*! \brief Get events that a state contains
     * \warning On large binery trees, it can be very expensive.
     *
     * @param aQ A state on the sys
     * \return Bit set with events of state events
     */
    EventsSet<NEvents> getStateEvents_impl(StorageIndex const& aQ) const;

    /*! \brief Get inverse events that a state contains
     *
     * @param aQ A state on the sys
     * \return Bit set with events index set to true
     */
    EventsSet<NEvents> getInvStateEvents_impl(StorageIndex const& aQ) const;

    /*! \brief Invert graph
     * \details This is used on some operations... it can be very inneficient
     * for very large graphs
     * It is const, since it changes only a mutable member
     *
     * \return void
     */
    void allocateInvertedGraph_impl() const;

    /*! \brief Free inverted graph
     * \details It is const, since it changes only a mutable member
     */
    void clearInvertedGraph_impl() const;

protected:
    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::SynchronizeStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::SynchronizeEmptyStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

    friend class cldes::op::SuperProxy<NEvents, StorageIndex>;

    /*! \brief Disabled default constructor
     * \details There is no use for the default constructor.
     */
    SyncSysProxy() = default;

private:
    friend DESystem SupervisorSynth<>(DESystemBase const& aP,
                                      DESystemBase const& aE,
                                      EventsTableHost const& aNonContr);

    /*! \brief Reference to the left operand
     */
    DESystemBase const& sys0_;

    /*! \brief Reference to right operand
     */
    DESystemBase const& sys1_;

    /*! \brief Cache number of states of Sys0
     *
     */
    StorageIndex n_states_sys0_;

    /*! \brief Virtual states contained in the current system
     */
    StatesTable virtual_states_;

    /*! \brief Events contained only in the left operator of a synchronizing op.
     */
    EventsSet<NEvents> only_in_0_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     */
    EventsSet<NEvents> only_in_1_;

    /*! \brief Pointer to real systems, if exists
     * \details Since it is a lazy system, it needs to be declared mutable for
     * enabling lazy evaluation.
     */
    std::shared_ptr<DESystem> sys_ptr_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     */
    TrVector transtriplet_;

    /*! \brief 3-tuples for filling graph_
     */
    std::vector<Triplet<NEvents>> triplet_;

    /*! \brief 3-tuples for filling bit_graph_
     */
    std::vector<BitTriplet> bittriplet_;
};

} // namespace op
} // namespace cldes

// include methods definitions
#include "cldes/src/operations/SyncSysProxyCore.hpp"

#endif // TRANSITION_PROXY_HPP
