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

 File: cldes/DESystemBase.hpp
 Description: DESystemBase abstract class definition.
 =========================================================================
*/

#ifndef DESYSTEMBASE_HPP
#define DESYSTEMBASE_HPP

#include "cldes/Constants.hpp"
#include "cldes/EventsSet.hpp"
#include <Eigen/Sparse>
#include <set>
#include <vector>

namespace cldes {

/*! \brief Vector of states type
 */
template<typename StorageIndex>
using StatesArray = std::vector<StorageIndex>;

/*! \brief Discrete-Events System Base
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystemBase
{
public:
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Set of states type
     */

    /*! \brief Set of states type
     */
    using StatesSet = std::set<StorageIndex>;

    /*! \brief Vector of states type
     */
    using StatesTable = std::vector<StorageIndex>;

    /*! \brief Table of transitions on a STL container
     */
    using StatesEventsTable = std::vector<EventsSet<NEvents>>;

    /*! \brief DESystem constructor with empty matrix
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     * @param aMarkedStates System's marked states
     */
    explicit inline DESystemBase(StorageIndex const& aStatesNumber,
                                 StorageIndex const& aInitState)
    {
        states_number_ = aStatesNumber;
        init_state_ = aInitState;
    };

    /*! \brief DESystem destructor
     */
    virtual ~DESystemBase() = default;

    /*! \brief Move constructor
     *
     * Enable move semantics
     */
    DESystemBase(DESystemBase&&) = default;

    /*! \brief Copy constructor
     *
     * Needs to define this, since move semantics is enabled
     */
    DESystemBase(DESystemBase const&) = default;

    /*! \brief Operator =
     *
     * Uses move semantics
     */
    DESystemBase& operator=(DESystemBase&&) = default;

    /*! \brief Operator = to const type
     *
     * Needs to define this, since move semantics is enabled
     */
    DESystemBase& operator=(DESystemBase const&) = default;

    /*! \brief Returns number of states of the system
     *
     * Returns states_value_ by value.
     */
    inline StorageIndex Size() const { return states_number_; }

    /*! \brief Returns events contained in the system
     *
     */
    inline EventsSet<NEvents> GetEvents() const { return events_; }

    /*! \brief Returns number of states contained in the system
     *
     */
    inline StorageIndex GetStatesNumber() const { return states_number_; }

    /*! \brief Returns number of states contained in the system
     *
     */
    inline StorageIndex GetInitialState() const { return init_state_; }

    /*! \brief Returns marked states
     *
     */
    inline StatesSet GetMarkedStates() const { return marked_states_; }

    /*! \brief Set system's number of states
     *
     */
    inline void SetStatesNumber(StorageIndex const& aStNum)
    {
        states_number_ = aStNum;
    }

    /*! \brief Set system's system events table
     *
     */
    inline void ResizeStatesEvents(StorageIndex const& aSize)
    {
        states_events_.resize(aSize);
        inv_states_events_.resize(aSize);
    }

    /*! \brief Set system's system events table
     *
     */
    inline void SetStateEvents(StorageIndex const& aQ,
                               EventsSet<NEvents> const& aEvents)
    {
        states_events_[aQ] = aEvents;
    }

    /*! \brief Set inverted states events
     *
     */
    inline void SetInvStateEvents(StorageIndex const& aQ,
                                  EventsSet<NEvents> const& aEvents)
    {
        inv_states_events_[aQ] = aEvents;
    }
    /*! \brief Set system's system events table
     *
     */
    inline void SetStatesEvents(StatesEventsTable const& aEvents)
    {
        states_events_ = aEvents;
    }

    /*! \brief Set inverted states events
     *
     */
    inline void SetInvStatesEvents(StatesEventsTable const& aEvents)
    {
        inv_states_events_ = aEvents;
    }

    /*! \brief Set inverted states events
     *
     */
    inline void SetEvents(EventsSet<NEvents> const& aEvents)
    {
        events_ = aEvents;
    }

    /*! \brief Set system's initial state
     *
     */
    inline void SetInitialState(StorageIndex const& aInitState)
    {
        init_state_ = aInitState;
    }

    /*! \brief Returns marked states
     *
     */
    inline void InsertMarkedState(StorageIndex const& aSt)
    {
        marked_states_.emplace(aSt);
    }

    /*! \brief Returns marked states
     *
     */
    inline void SetMarkedStates(StatesSet const& aStSet)
    {
        marked_states_ = aStSet;
    }

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    virtual bool ContainsTrans(StorageIndex const& aQ,
                               ScalarType const& aEvent) const = 0;

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    virtual StorageIndexSigned Trans(StorageIndex const& aQ,
                                     ScalarType const& aEvent) const = 0;

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    virtual bool ContainsInvTrans(StorageIndex const& aQ,
                                  ScalarType const& aEvent) const = 0;

    /*! \brief Returns DES inverse transition: q = f^-1(q_to, e)
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    virtual StatesArray<StorageIndex> InvTrans(
      StorageIndex const& aQfrom,
      ScalarType const& aEvent) const = 0;
    /*! \brief Returns EventsSet relative to state q
     *
     * @param aQ A state on the sys
     */
    virtual EventsSet<NEvents> GetStateEvents(StorageIndex const& aQ) const = 0;

    /*! \brief Returns EventsSet relative to state inv q
     *
     * @param aQ A state on the sys
     */
    virtual EventsSet<NEvents> GetInvStateEvents(
      StorageIndex const& aQ) const = 0;

    /*! \brief Invert graph
     *
     */
    virtual void AllocateInvertedGraph() const = 0;

    /*! \brief Free inverted graph
     *
     */
    virtual void ClearInvertedGraph() const = 0;

protected:
    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    DESystemBase() = default;

    /*! \brief Current system's states number
     *
     * Hold the number of states that the automata contains. As the automata
     * can be cut, the states number is not a constant at all.
     */
    StorageIndex states_number_;

    /*! \brief Current system's initial state
     *
     * Hold the initial state position.
     */
    StorageIndex init_state_;

    /*! \brief System's events hash table
     *
     * A hash table containing all the events that matter for the current
     * system.
     */
    EventsSet<NEvents> events_;

    /*! \brief Current system's marked states
     *
     * Hold all marked states. Cannot be const, since the automata can be
     * cut, and some marked states may be deleted.
     */
    StatesSet marked_states_;

    /*! \brief Vector containing a events hash table per state
     */
    StatesEventsTable states_events_;

    /*! \brief Vector containing a events hash table per state
     *
     * It represents the transitions of the inverted graph for the supervisor
     * synthesis.
     */
    StatesEventsTable inv_states_events_;
};
}
#endif // DESYSTEMBASE_HPP
