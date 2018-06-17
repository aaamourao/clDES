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
#include "cldes/src/operations/SyncSysProxyCore.hpp"

namespace cldes {
namespace op {

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeEmptyStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

// Alias to events hash map
using EventsTableHost = spp::sparse_hash_set<uint8_t>;

// Forward declaration of friend function
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(DESystemBase<NEvents, StorageIndex> const& aP,
                DESystemBase<NEvents, StorageIndex> const& aE,
                EventsTableHost const& aNonContr);

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
class SyncSysProxy : public DESystemBase<NEvents, StorageIndex>
{
public:
    /*! \brief Signed template parameter type for eigen indexes
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Base alias
     * \details Alias to implicit speciallization of base class.
     */
    using DESystemBase = DESystemBase<NEvents, StorageIndex>;

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
    inline SyncSysProxy(DESystemBase const& aSys0, DESystemBase const& aSys1)
      : DESystemBase{ aSys0.GetStatesNumber() * aSys1.GetStatesNumber(),
                      aSys1.GetInitialState() * aSys0.GetStatesNumber() +
                        aSys0.GetInitialState() }
      , sys0_{ aSys0 }
      , sys1_{ aSys1 }
    {
        sys_ptr_ = nullptr;
        n_states_sys0_ = aSys0.GetStatesNumber();

        auto const in_both = aSys0.GetEvents() & aSys1.GetEvents();

        only_in_0_ = aSys0.GetEvents() ^ in_both;
        only_in_1_ = aSys1.GetEvents() ^ in_both;

        // Initialize events_ from base class
        this->events_ = aSys0.GetEvents() | aSys1.GetEvents();

        for (auto q0 : aSys0.GetMarkedStates()) {
            for (auto q1 : aSys1.GetMarkedStates()) {
                this->marked_states_.emplace(q1 * n_states_sys0_ + q0);
            }
        }
    }

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
    operator DESystem()
    {
        if (virtual_states_.empty()) {
            SynchronizeEmptyStage2(*this);
        } else {
            std::sort(virtual_states_.begin(), virtual_states_.end());
            SynchronizeStage2(*this);
        }

        // Allocate memory for the real sys
        sys_ptr_ = std::make_shared<DESystem>(DESystem{});

        sys_ptr_->states_number_ = this->states_number_;
        sys_ptr_->init_state_ = this->init_state_;
        sys_ptr_->marked_states_ = this->marked_states_;
        sys_ptr_->states_events_ = this->states_events_;
        sys_ptr_->inv_states_events_ = this->inv_states_events_;
        sys_ptr_->events_ = this->events_;

        // Resize adj matrices
        sys_ptr_->graph_.resize(this->states_number_, this->states_number_);
        sys_ptr_->bit_graph_.resize(this->states_number_, this->states_number_);

        // Move triplets to graph storage
        sys_ptr_->graph_.setFromTriplets(triplet_.begin(), triplet_.end());
        sys_ptr_->bit_graph_.setFromTriplets(
          bittriplet_.begin(), bittriplet_.end(), [](bool const&, bool const&) {
              return true;
          });

        triplet_.clear();
        bittriplet_.clear();

        sys_ptr_->graph_.makeCompressed();
        sys_ptr_->bit_graph_.makeCompressed();

        return *sys_ptr_;
    }

    /*! \brief Convert to DESystem from const proxy
     * \details Convert it to a concrete system when the current system is
     * const.
     * \warning Expensive, because it implies a copy. Avoid it.
     */
    operator DESystem() const
    {
        auto cp = *this;
        return DESystem(cp);
    }

    /*! \brief Is it real?
     * \details Nooo!!!
     *
     * \return False, always.
     */
    inline bool IsVirtual() const override
    {
        return true;
    }

    /*! \brief Clone method for polymorphic copy
     * \details Method for cloning on a polymorphic way.
     *
     * \return shared pointer to base.
     */
    inline std::shared_ptr<DESystemBase> Clone() const override
    {
        std::shared_ptr<DESystemBase> this_ptr =
          std::make_shared<SyncSysProxy>(*this);
        return this_ptr;
    }

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    inline bool ContainsTrans(StorageIndex const& aQ,
                              ScalarType const& aEvent) const override
    {
        // q = (qx, qy)
        auto const qx = aQ % n_states_sys0_;
        auto const qy = aQ / n_states_sys0_;

        auto const in_x = sys0_.ContainsTrans(qx, aEvent);
        auto const in_y = sys1_.ContainsTrans(qy, aEvent);

        auto contains = false;

        if ((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
            (in_y && only_in_1_.test(aEvent))) {
            contains = true;
        }

        return contains;
    }

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    inline StorageIndexSigned Trans(StorageIndex const& aQ,
                                    ScalarType const& aEvent) const override
    {
        // q = (qx, qy)
        auto const qx = aQ % n_states_sys0_;
        auto const qy = aQ / n_states_sys0_;

        auto const in_x = sys0_.ContainsTrans(qx, aEvent);
        auto const in_y = sys1_.ContainsTrans(qy, aEvent);

        if (!((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
              (in_y && only_in_1_.test(aEvent)))) {
            return -1;
        }

        if (in_x && in_y) {
            auto const q0 = sys0_.Trans(qx, aEvent);
            auto const q1 = sys1_.Trans(qy, aEvent);

            return q1 * n_states_sys0_ + q0;
        } else if (in_x) {
            auto const q0 = sys0_.Trans(qx, aEvent);
            return qy * n_states_sys0_ + q0;
        } else { // in_y
            auto const q1 = sys1_.Trans(qy, aEvent);
            return q1 * n_states_sys0_ + qx;
        }
    }

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    inline bool ContainsInvTrans(StorageIndex const& aQ,
                                 ScalarType const& aEvent) const override
    {
        // q = (qx, qy)
        auto const qx = aQ % n_states_sys0_;
        auto const qy = aQ / n_states_sys0_;

        auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
        auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

        auto contains = false;

        if ((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
            (in_y && only_in_1_.test(aEvent))) {
            contains = true;
        }

        return contains;
    }

    /*! \brief Returns DES inverse transition
     * \details  q = f^-1(q_to, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    inline StatesArray<StorageIndex> InvTrans(
      StorageIndex const& aQ,
      ScalarType const& aEvent) const override
    {
        // q = (qx, qy)
        auto const qx = aQ % n_states_sys0_;
        auto const qy = aQ / n_states_sys0_;

        auto const in_x = sys0_.ContainsInvTrans(qx, aEvent);
        auto const in_y = sys1_.ContainsInvTrans(qy, aEvent);

        StatesArray<StorageIndex> inv_transitions;

        if (!((in_x && in_y) || (in_x && only_in_0_.test(aEvent)) ||
              (in_y && only_in_1_.test(aEvent)))) {
            return inv_transitions;
        }

        if (in_x && in_y) {
            auto const inv_trans_0 = sys0_.InvTrans(qx, aEvent);
            auto const inv_trans_1 = sys1_.InvTrans(qy, aEvent);

            inv_transitions.reserve(inv_trans_0.size() + inv_trans_1.size());

            for (auto q0 : inv_trans_0) {
                for (auto q1 : inv_trans_1) {
                    auto const q_from = q1 * n_states_sys0_ + q0;
                    inv_transitions.push_back(q_from);
                }
            }
        } else if (in_x) {
            auto const inv_trans_0 = sys0_.InvTrans(qx, aEvent);

            inv_transitions.reserve(inv_trans_0.size());

            for (auto q : inv_trans_0) {
                auto const q_from = qy * n_states_sys0_ + q;
                inv_transitions.push_back(q_from);
            }
        } else { // in_y
            auto const inv_trans_1 = sys1_.InvTrans(qy, aEvent);

            inv_transitions.reserve(inv_trans_1.size());

            for (auto q : inv_trans_1) {
                auto const q_from = q * n_states_sys0_ + qx;
                inv_transitions.push_back(q_from);
            }
        }

        return inv_transitions;
    }

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
    inline EventsSet<NEvents> GetStateEvents(
      StorageIndex const& aQ) const override
    {
        auto const state_event_0 = sys0_.GetStateEvents(aQ % n_states_sys0_);
        auto const state_event_1 = sys1_.GetStateEvents(aQ / n_states_sys0_);
        auto const state_event = (state_event_0 & state_event_1) |
                                 (state_event_0 & only_in_0_) |
                                 (state_event_1 & only_in_1_);

        return state_event;
    }

    /*! \brief Get inverse events that a state contains
     *
     * @param aQ A state on the sys
     * \return Bit set with events index set to true
     */
    inline EventsSet<NEvents> GetInvStateEvents(
      StorageIndex const& aQ) const override
    {
        auto const state_event_0 = sys0_.GetInvStateEvents(aQ % n_states_sys0_);
        auto const state_event_1 = sys1_.GetInvStateEvents(aQ / n_states_sys0_);
        auto const state_event = (state_event_0 & state_event_1) |
                                 (state_event_0 & only_in_0_) |
                                 (state_event_1 & only_in_1_);

        return state_event;
    }

    /*! \brief Invert graph
     * \details This is used on some operations... it can be very inneficient
     * for very large graphs
     * It is const, since it changes only a mutable member
     *
     * \return void
     */
    inline void AllocateInvertedGraph() const override
    {
        sys0_.AllocateInvertedGraph();
        sys1_.AllocateInvertedGraph();
    }

    /*! \brief Free inverted graph
     * \details It is const, since it changes only a mutable member
     */
    inline void ClearInvertedGraph() const override
    {
        sys0_.ClearInvertedGraph();
        sys1_.ClearInvertedGraph();
    }

protected:
    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::SynchronizeStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

    /*! \brief Second step of the lazy parallel composition
     */
    friend void cldes::op::SynchronizeEmptyStage2<>(
      SyncSysProxy<NEvents, StorageIndex>& aVirtualSys);

    /*! \brief Monolithic supervisor synthesis
     */
    friend DESystem SupervisorSynth<>(DESystemBase const& aP,
                                      DESystemBase const& aE,
                                      EventsTableHost const& aNonContr);

    /*! \brief Disabled default constructor
     * \details There is no use for the default constructor.
     */
    inline SyncSysProxy() = default;

private:
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

#endif // TRANSITION_PROXY_HPP