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
    using Base =
      DESystemBase<NEvents, StorageIndex, SyncSysProxy<NEvents, StorageIndex>>;

    /*! \brief Alias to the the related DESystem template
     * \details alias Alias to implicit speciallization to concrete system.
     */
    using RealSys = DESystem<NEvents, StorageIndex>;
    using BaseReal = DESystemBase<NEvents, StorageIndex, RealSys>;

    /*! \brief Vector of states type
     * \details Vector containing states represented by unsigned indexes.
     */
    using StatesTable = typename BaseReal::StatesTable;

    /*! \brief Vector of inverted transitions
     * \details Vector which stores transitions:
     * f(s, e) = s_out -> (s_out, (s, e)) is the inverted transition.
     */
    using TrVector =
      std::vector<std::pair<StorageIndex, InvArgtrans<StorageIndex>*>>;

    /*! \brief SyncSysProxy unique constructor
     * Create a binary tree that represents multiple operations.
     *
     * @param aSys0 Left operand DESystem reference
     * @param aSys1 Right operand DESystem reference
     */
    SyncSysProxy(BaseReal const& aSys0, BaseReal const& aSys1);

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
    operator RealSys() noexcept;

    /*! \brief Is it real?
     * \details Nooo!!!
     *
     * \return False, always.
     */
    bool constexpr static isVirtual_impl() noexcept { return true; }

    /*! \brief clone method for polymorphic copy
     * \details Method for cloning on a polymorphic way.
     *
     * \return shared pointer to base.
     */
    std::shared_ptr<Base> constexpr clone_impl() const noexcept
    {
        std::shared_ptr<Base> this_ptr = std::make_shared<SuperProxy>(*this);
        return this_ptr;
    }

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool containstrans_impl(StorageIndex const& aQ,
                            ScalarType const& aEvent) const noexcept;

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StorageIndexSigned trans_impl(StorageIndex const& aQ,
                                  ScalarType const& aEvent) const noexcept;

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool containsinvtrans_impl(StorageIndex const& aQ,
                               ScalarType const& aEvent) const;

    /*! \brief Returns DES inverse transition
     * \details  q = f^-1(q_to, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StatesArray<StorageIndex> invtrans_impl(StorageIndex const& aQ,
                                            ScalarType const& aEvent) const;

    /*! \brief Get events that a state contains
     * \warning On large binery trees, it can be very expensive.
     *
     * @param aQ A state on the sys
     * \return Bit set with events of state events
     */
    EventsSet<NEvents> getStateEvents_impl(StorageIndex const& aQ) const
      noexcept;

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
    void allocateInvertedGraph_impl() const noexcept;

    /*! \brief Free inverted graph
     * \details It is const, since it changes only a mutable member
     */
    void clearInvertedGraph_impl() const noexcept;

    /*! \brief Observer property checker
     *
     * @param aLang Language
     */
    bool constexpr checkObsProp_impl(EventsSet<NEvents> const&) const noexcept
    {
        return false;
    }

    /*! \brief Projection operation
     *
     * @param aAlphanet Events set which the system will be projected.
     */
    SyncSysProxy constexpr& proj_impl(EventsSet<NEvents> const&) noexcept
    {
        return *this;
    }

    /*! \brief Inverse projection operation
     *
     * @param aLang Set of events
     */
    SyncSysProxy constexpr& invproj_impl(EventsSet<NEvents> const&) const
      noexcept
    {
        return *this;
    }

protected:
    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::synchronizeStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys) noexcept;

    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::synchronizeEmptyStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys) noexcept;

    friend class cldes::op::SuperProxy<NEvents, StorageIndex>;

    /*! \brief Disabled default constructor
     * \details There is no use for the default constructor.
     */
    SyncSysProxy() = default;

private:
    friend RealSys supC<>(BaseReal const& aP,
                          BaseReal const& aE,
                          EventsTableHost const& aNonContr) noexcept;

#ifdef __GNUC__
    friend transMap<StorageIndex> computeSupCStates_<>(
      SyncSysProxy<NEvents, StorageIndex> const& aVirtualSys,
      EventsSet<NEvents> const&& aNonContrBit,
      EventsSet<NEvents> const&& aPNonContrBit,
      BaseReal const& aP) noexcept;

    friend void processVirtSys_<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys,
      unsigned long const& aSparcityPattern,
      SparseStatesMap<StorageIndex>&& aStatesMap) noexcept;
#elif __clang__
    friend inline transMap<StorageIndex> computeSupCStates_<>(
      SyncSysProxy<NEvents, StorageIndex> const& aVirtualSys,
      EventsSet<NEvents> const&& aNonContrBit,
      EventsSet<NEvents> const&& aPNonContrBit,
      BaseReal const& aP) noexcept;

    friend inline void processVirtSys_<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys,
      unsigned long const& aSparcityPattern,
      SparseStatesMap<StorageIndex>&& aStatesMap) noexcept;
#endif

    friend StorageIndex aproxSpacPat_<>(
      SyncSysProxy<NEvents, StorageIndex> const& aV) noexcept;

    /*! \brief Reference to the left operand
     */
    BaseReal const& sys0_;

    /*! \brief Reference to right operand
     */
    BaseReal const& sys1_;

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

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     */
    TrVector transtriplet_;

    /*! \brief 3-tuples for filling graph_
     */
    std::vector<Triplet<NEvents>> triplet_;
};

} // namespace op
} // namespace cldes

// include methods definitions
#include "cldes/src/operations/SyncSysProxyCore.hpp"

#endif // TRANSITION_PROXY_HPP
