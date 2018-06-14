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

#ifndef SYNC_SYS_PROXY_HPP
#define SYNC_SYS_PROXY_HPP

#include "cldes/Constants.hpp"
#include "cldes/DESystemBase.hpp"
#include "cldes/EventsSet.hpp"

namespace cldes {
namespace op {

/*! \brief Proxy to a virtual sync sys
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type for storing the indexes of each state
 */
template<uint8_t NEvents, typename StorageIndex>
class SyncSysProxy : public DESystemBase<NEvents, StorageIndex>
{
public:
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Vector of states type
     */
    using StatesTable =
      typename DESystemBase<NEvents, StorageIndex>::StatesTable;

    /*! \brief Vector of inverted transitions
     *
     * f(s, e) = s_out -> (s_out, (s, e)) is the inverted transition.
     */
    using TrVector =
      std::vector<std::pair<StorageIndex, InvArgTrans<StorageIndex>*>>;

    /*! \brief SyncSysProxy unique constructor
     *
     * Feed const data-members.
     *
     * @param aSys0Ptr Left operand DESystem reference
     * @param aSys1Ptr Right operand DESystem reference
     */
    explicit SyncSysProxy(DESystemBase<NEvents, StorageIndex> const& aSys0,
                          DESystemBase<NEvents, StorageIndex> const& aSys1);

    /*! \brief DESystem destructor
     */
    virtual ~SyncSysProxy() = default;

    /*! \brief Move constructor
     *
     * Enable move semantics
     */
    SyncSysProxy(SyncSysProxy&&) = default;

    /*! \brief Copy constructor
     *
     * Needs to define this, since move semantics is enabled
     */
    SyncSysProxy(SyncSysProxy const&) = default;

    /*! \brief Operator =
     *
     * Uses move semantics
     */
    SyncSysProxy<NEvents, StorageIndex>& operator=(SyncSysProxy&&) = default;

    /*! \brief Operator = to const type
     *
     * Needs to define this, since move semantics is enabled
     */
    SyncSysProxy<NEvents, StorageIndex>& operator=(SyncSysProxy const&) =
      default;

    /*! \brief Overload conversion to DESystem
     *
     * @param aQ State
     * @param aEvent Event
     */
    operator DESystem<NEvents, StorageIndex>();

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool ContainsTrans(StorageIndex const& aQ,
                       ScalarType const& aEvent) const override;

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StorageIndexSigned Trans(StorageIndex const& aQ,
                             ScalarType const& aEvent) override;

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    bool ContainsInvTrans(StorageIndex const& aQ,
                          ScalarType const& aEvent) const override;

    /*! \brief Returns DES inverse transition: q = f^-1(q_to, e)
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    StatesArray<StorageIndex> InvTrans(StorageIndex const& aQfrom,
                                       ScalarType const& aEvent) override;

protected:
    /*! \brief Disabled default constructor
     *
     * There is no use for the default constructor.
     */
    explicit SyncSysProxy();

private:
    /*! \brief Raw pointer to DESystemBase object
     *
     * Raw pointer to the owner of the proxied element.
     */
    DESystemBase<NEvents, StorageIndex>& sys0_;

    /*! \brief Raw pointer to DESystemBase object
     *
     * Raw pointer to the owner of the proxied element.
     */
    DESystemBase<NEvents, StorageIndex>& sys1_;

    /*! \brief Cache number of states of Sys0
     *
     */
    StorageIndex n_states_sys0_;

    /*! \brief Virtual states contained in the current system
     *
     */
    StatesTable virtual_states_;

    /*! \brief Events contained only in the left operator of a synchronizing op.
     *
     */
    EventsSet<NEvents> only_in_0_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     *
     */
    EventsSet<NEvents> only_in_1_;

    /*! \brief Events contained only in the right operator of a synchronizing
     * op.
     *
     */
    TrVector transtriplet_;
};

} // namespace op
} // namespace cldes

// including implementation
#include "cldes/src/operations/SyncSysProxyCore.hpp"

#endif // TRANSITION_PROXY_HPP
