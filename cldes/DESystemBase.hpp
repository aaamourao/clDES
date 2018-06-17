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
/*!
 * \file cldes/DESystemBase.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-06-16
 *
 * DESystemBase template abstract class declaration and definition .
 */

#ifndef DESYSTEMBASE_HPP
#define DESYSTEMBASE_HPP

#include "cldes/Constants.hpp"
#include "cldes/EventsSet.hpp"
#include "cldes/src/des/DESystemBaseCore.hpp"

namespace cldes {

/*! \brief Vector of states type
 * \details Vector containing states represented by usigned numbers.
 */
template<typename StorageIndex>
using StatesArray = std::vector<StorageIndex>;

/*! \class DESystemBase
 *  \brief A discrete-events system base abstract class
 *  \details DESystemBase implements discrete-events interface and
 *  its basic properties.
 *
 * \tparam NEvents Number of events
 * \tparam StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystemBase
{
public:
    /*! \brief Signed type of indexes
     * \details Eigen uses signed value for indexes.
     */
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;

    /*! \brief Vector of states type
     *  \details Vector of usigned interger type which represent states.
     */
    using StatesSet = std::set<StorageIndex>;

    /*! \brief Table of transitions on a STL container
     *  \details Vector containing bitsets which represent the events that a
     *  state represented by the vector index contains.
     */
    using StatesTable = std::vector<StorageIndex>;

    /*! \brief Table of transitions on a STL container
     *  \details Vector containing bitsets which represent the events that a
     *  state represented by the vector index contains.
     */
    using StatesEventsTable = std::vector<EventsSet<NEvents>>;

    /*! \brief DESystem constructor
     * \details Create a system with the defined params.
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     */
    inline DESystemBase(StorageIndex const& aStatesNumber,
                        StorageIndex const& aInitState)
    {
        states_number_ = aStatesNumber;
        init_state_ = aInitState;
    }

    /*! \brief Default constructor
     * \details Creates an empty system with 0 states and initial
     * states 0. Marked states is undefined.
     *
     */
    inline DESystemBase()
    {
        states_number_ = 0;
        init_state_ = 0;
        events_ = EventsSet<NEvents>{};
        marked_states_ = StatesSet{};
    };

    /*! \brief DESystem destructor
     * \details It is virtual for avoiding issues with polymorphism and memory
     * leaks.
     */
    virtual ~DESystemBase() = default;

    /*! \brief Move constructor
     * \details Enable move semantics.
     */
    DESystemBase(DESystemBase&&) = default;

    /*! \brief Copy constructor
     * \details Need to define this to enable copy by reference.
     */
    DESystemBase(DESystemBase const&) = default;

    /*! \brief Overload operator =
     * \details Use move semantics when assigning to rvalues.
     */
    DESystemBase& operator=(DESystemBase&&) = default;

    /*! \brief Copy constructor
     * \details Need to define this to enable copy by reference.
     */
    DESystemBase& operator=(DESystemBase const&) = default;

    /*! \brief Returns number of states of the system
     * \details Returns states_value_ by value.
     *
     * \return The number of states contained in the current system
     */
    inline StorageIndex Size() const { return states_number_; }

    /*! \brief Get number of states of the current system
     * \details states_number_ getter.
     *
     * \return The number of states contained in the current system
     */
    inline EventsSet<NEvents> GetEvents() const { return events_; }

    /*! \brief Returns number of states contained in the system
     *
     * \return Unsigned integer type represent the system's states number
     */
    inline StorageIndex GetStatesNumber() const { return states_number_; }

    /*! \brief Returns number of states contained in the system
     *
     * \return Unsigned integer type representing the initial state
     */
    inline StorageIndex GetInitialState() const { return init_state_; }

    /*! \brief Returns marked states
     *
     * \return Set of usigned integer type representing the marked states.
     */
    inline StatesSet GetMarkedStates() const { return marked_states_; }

    /*! \brief Set inverted states events
     *
     * \param aEvents Bit set with new events of the system
     * \return void
     */
    inline void SetEvents(EventsSet<NEvents> const& aEvents)
    {
        events_ = aEvents;
    }

    /*! \brief Set system's number of states
     *
     * \param aStNum New system's states number.
     * \return void
     */
    inline void SetStatesNumber(StorageIndex const& aStNum)
    {
        states_number_ = aStNum;
    }

    /*! \brief Set system's initial state
     *
     * \param aInitState New initial state
     * \return void
     */
    inline void SetInitialState(StorageIndex const& aInitState)
    {
        init_state_ = aInitState;
    }

    /*! \brief Returns marked states
     *
     * \param aSt state which will be inserted
     * \return void
     */
    void InsertMarkedState(StorageIndex const& aSt)
    {
        marked_states_.emplace(aSt);
    }

    /*! \brief Set system's marked states
     *
     * \param aStSet Set of states
     * \return void
     */
    inline void SetMarkedStates(StatesSet const& aStSet)
    {
        marked_states_ = aStSet;
    }

    /*! \brief Resize state_events
     * \details Necessary to run it when inserting or remove
     * events.
     * \warning It is used by many operations. Let them do it for you.
     * You can corrupt a system by executing it.
     *
     * \param aSize New state_events_ size
     * \return void
     */
    inline void ResizeStatesEvents(StorageIndex const& aSize)
    {
        states_events_.resize(aSize);
        inv_states_events_.resize(aSize);
    }

    /*! \brief Resize state_events
     * \details Necessary to run it when inserting or remove
     * events.
     * \warning It is used by many operations. Let them do it for you.
     * You can corrupt a system by executing it.
     *
     * \param aEvents Events set vector
     * \return void
     */
    inline void SetStatesEvents(StatesEventsTable const& aEvents)
    {
        states_events_ = aEvents;
    }

    /*! \brief Set inv_state_events
     * \details It is used by many operations. Set the inverse events
     * of each states contains.
     * \warning It is used by many operations. Let them do it for you
     * You can corrupt a system by executing it.
     *
     * \param aEvents Events set vector
     * \return void
     */
    inline void SetInvStatesEvents(StatesEventsTable const& aEvents)
    {
        inv_states_events_ = aEvents;
    }

    /*! \brief Set state_events of a specific state
     * \details It is used by many operations. Set the events of
     * a state contains.
     * \warning It is used by many operations. Let them do it for you
     * You can corrupt a system by executing it.
     *
     * \param aQ State represented by an unsigned integer type
     * \param aEvent A event represented by a 8 bit unsigned integer
     * \return void
     */
    inline void SetStateEvents(StorageIndex const& aQ,
                               EventsSet<NEvents> const& aEvent)
    {
        states_events_[aQ] = aEvent;
    }

    /*! \brief Set inv_state_events of a specific state
     * \details It is used by many operations. Set the events of
     * a state contains.
     * \warning It is used by many operations. Let them do it for you
     * You can corrupt a system by executing it.
     *
     * \param aQ State represented by an unsigned integer type
     * \param aEvent A event represented by a 8 bit unsigned integer
     * \return void
     */
    inline void SetInvStateEvents(StorageIndex const& aQ,
                                  EventsSet<NEvents> const& aEvent)
    {
        inv_states_events_[aQ] = aEvent;
    }

    /*! \brief Is it real?
     *
     * \return Boolean with the answer. It is true or false, not 42.
     */
    virtual bool IsVirtual() const = 0;

    /*! \brief Clone method to enable poliphormism
     *
     * \return shared pointer of type base to the system
     */
    virtual std::shared_ptr<DESystemBase> Clone() const = 0;

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
     * @param aQ State
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
