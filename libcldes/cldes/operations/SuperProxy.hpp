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

 File: cldes/operations/SuperSynth.hpp
 Description: Virtual Proxy for an implicit monolithic supervisor
 synthesis
 =========================================================================
*/
/*!
 * \file cldes/operations/SuperSynth.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-07-11
 *
 * Virtual Proxy for an implicit monolithic supervisor synthesis.
 */

#ifndef SUPER_PROXY_HPP
#define SUPER_PROXY_HPP

#include "cldes/Constants.hpp"
#include "cldes/DESystemBase.hpp"
#include "cldes/EventsSet.hpp"
#include "cldes/src/operations/SuperProxyFwd.hpp"

namespace cldes {
namespace op {

/*! \class SuperProxy
 * \brief Monolothic supervisor synthesis return type
 * \details Represents a supervisor synthesis. Allows lazy evaluation, but can
 * be converted to a concrete system any time with a DESystem() cast.
 *
 * \tparam NEvents Number of events
 * \tparam StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<class SysT_l, class SysT_r>
class SuperProxy
  : public DESystemBase<SysTraits<SysT_l>::Ne_,
                        typename SysTraits<SysT_l>::Si_,
                        SuperProxy<SysT_l, SysT_r>>
{
public:
    uint8_t static constexpr NEvents = SysTraits<SysT_l>::Ne_;
    using StorageIndex = typename SysTraits<SysT_l>::Si_;

    using EventsSet_t = EventsSet<NEvents>;
    /*! \brief Signed template parameter type for eigen indexes
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Base alias
     * \details Alias to implicit speciallization of base class.
     */
    using Base =
      DESystemBase<NEvents, StorageIndex, SuperProxy<SysT_l, SysT_r>>;

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
      std::vector<std::pair<StorageIndex, InvArgtrans<StorageIndex>>>;

    /*! \brief SyncSysProxy unique constructor
     * Create a binary tree that represents multiple operations.
     *
     * @param aPlant Left operand DESystem reference
     * @param aSpec Right operand DESystem reference
     */
    SuperProxy(SysT_l const& aPlant,
               SysT_r const& aSpec,
               EventsTableHost const& aNonContr);

    /*! \brief DESystem destructor
     * \details Override base destructor.
     */
    ~SuperProxy() = default;

    /*! \brief Move constructor
     * \details Enable move semantics.
     */
    SuperProxy(SuperProxy&&) = default;

    /*! \brief Copy constructor
     * \details Enable copy by reference.
     */
    SuperProxy(SuperProxy const&) = default;

    /*! \brief Operator =
     * Use move semantics when assigning to rvalues.
     */
    SuperProxy& operator=(SuperProxy&&) = default;

    /*! \brief Operator = to const type
     * \details Enables copy by reference.
     */
    SuperProxy& operator=(SuperProxy const&) = default;

    /*! \brief Convert to DESystem from const proxy
     * \details Convert it to a concrete system when the current system is
     * const.
     * \warning Expensive, because it implies a copy. Avoid it.
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

    void trim() noexcept;

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
    SuperProxy constexpr& proj_impl(EventsSet<NEvents> const&) noexcept
    {
        return *this;
    }

    /*! \brief Inverse projection operation
     *
     * @param aLang Set of events
     */
    SuperProxy constexpr& invproj_impl(EventsSet<NEvents> const&) const noexcept
    {
        return *this;
    }

protected:
    /*! \brief Monolithic supervisor synthesis
     */
    friend RealSys supC<>(SysT_l const& aP,
                          SysT_r const& aE,
                          EventsTableHost const& aNonContr) noexcept;

    /*! \brief Disabled default constructor
     * \details There is no use for the default constructor.
     */
    SuperProxy() = default;

    /*! \brief Find states that are not in SupC
     */
    void findRemovedStates_(SysT_l const& aP,
                            SysT_r const& aE,
                            EventsTableHost const& aNonContr) noexcept;

private:
    /*! \brief Reference to the left operand
     */
    SysT_l const& sys0_;

    /*! \brief Reference to right operand
     */
    SysT_r const& sys1_;

    /*! \brief Cache number of states of Sys0
     *
     */
    StorageIndex n_states_sys0_;

    /*! \brief Virtual states contained in the current system
     */
    StatesTableHost<StorageIndex> virtual_states_;

    /*! \brief Events contained only in the left operator of a synchronizing op.
     */
    EventsSet<NEvents> only_in_plant_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     */
    EventsSet<NEvents> only_in_spec_;

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
#include "cldes/src/operations/SuperProxyCore.hpp"

#endif // SUPER_PROXY_HPP
