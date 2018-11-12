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
#include "cldes/src/des/DESystemBaseFwd.hpp"

namespace cldes {

/*! \class DESystemBase
 *  \brief A discrete-events system base abstract class
 *  \details DESystemBase implements discrete-events interface and
 *  its basic properties. It is implemented using CRTP pattern, static
 *  polymorphism, for the sake of speed.
 *
 * \tparam NEvents Number of events
 * \tparam StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
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
    StorageIndex constexpr size() const noexcept { return states_number_; }

    /*! \brief Get number of states of the current system
     * \details states_number_ getter.
     *
     * \return The number of states contained in the current system
     */
    EventsSet<NEvents> constexpr getEvents() const noexcept { return events_; }

    /*! \brief Returns number of states contained in the system
     *
     * \return Unsigned integer type represent the system's states number
     */
    StorageIndex constexpr getStatesNumber() const noexcept
    {
        return states_number_;
    }

    /*! \brief Returns number of states contained in the system
     *
     * \return Unsigned integer type representing the initial state
     */
    StorageIndex constexpr getInitialState() const noexcept
    {
        return init_state_;
    }

    /*! \brief Returns marked states
     *
     * \return Set of usigned integer type representing the marked states.
     */
    StatesSet constexpr getMarkedStates() const noexcept
    {
        return marked_states_;
    }

    /*! \brief Set inverted states events
     *
     * \param aEvents Bit set with new events of the system
     * \return void
     */
    void setEvents(EventsSet<NEvents> const& aEvents) noexcept;

    /*! \brief Set system's number of states
     *
     * \param aStNum New system's states number.
     * \return void
     */
    void setStatesNumber(StorageIndex const& aStNum) noexcept;

    /*! \brief Set system's initial state
     *
     * \param aInitState New initial state
     * \return void
     */
    void setInitialState(StorageIndex const& aInitState) noexcept;

    /*! \brief Returns marked states
     *
     * \param aSt state which will be inserted
     * \return void
     */
    void insertMarkedState(StorageIndex const& aSt) noexcept;

    /*! \brief Set system's marked states
     *
     * \param aStSet Set of states
     * \return void
     */
    void setMarkedStates(StatesSet const& aStSet) noexcept;

    /*! \brief Resize state_events
     * \details Necessary to run it when inserting or remove
     * events.
     * \warning It is used by many operations. Let them do it for you.
     * You can corrupt a system by executing it.
     *
     * \param asize New state_events_ size
     * \return void
     */
    void resizeStatesEvents(StorageIndex const& asize) noexcept;

    /*! \brief Force statesevents to assume value
     * \details Necessary to run it when inserting or remove
     * events.
     * \warning It is used by many operations. Let them do it for you.
     * You can corrupt a system by executing it.
     *
     * \param aEvents Events set vector
     * \return void
     */
    void setStatesEvents(StatesEventsTable const& aEvents) noexcept;

    /*! \brief Set inv_state_events
     * \details It is used by many operations. Set the inverse events
     * of each states contains.
     * \warning It is used by many operations. Let them do it for you
     * You can corrupt a system by executing it.
     *
     * \param aEvents Events set vector
     * \return void
     */
    void setInvStatesEvents(StatesEventsTable const& aEvents) noexcept;

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
    void setStateEvents(StorageIndex const& aQ,
                        EventsSet<NEvents> const& aEvent) noexcept;

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
    void setInvStateEvents(StorageIndex const& aQ,
                           EventsSet<NEvents> const& aEvent) noexcept;

    /*! \brief Is it real?
     *
     * \return Boolean with the answer. It is true or false, not 42.
     */
    bool constexpr isVirtual() const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.isVirtual_impl();
    }

    /*! \brief clone method to enable poliphormism
     *
     * \return shared pointer of type base to the system
     */
    std::shared_ptr<DESystemBase> constexpr clone() const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.clone_impl();
    }

    /*! \brief Returns true if DES transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool constexpr containstrans(StorageIndex const& aQ,
                                 ScalarType const& aEvent) const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.containstrans_impl(aQ, aEvent);
    }

    /*! \brief Returns DES transition: q_to = f(q, e)
     *
     * @param aQ State
     * @param aEvent Event
     */
    StorageIndexSigned constexpr trans(StorageIndex const& aQ,
                                       ScalarType const& aEvent) const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.trans_impl(aQ, aEvent);
    }

    /*! \brief Returns true if DES inverse transition exists
     *
     * @param aQ State
     * @param aEvent Event
     */
    bool constexpr containsinvtrans(StorageIndex const& aQ,
                                    ScalarType const& aEvent) const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.containsinvtrans_impl(aQ, aEvent);
    }

    /*! \brief Returns DES inverse transition: q = f^-1(q_to, e)
     *
     * @param aQfrom State
     * @param aEvent Event
     */
    StatesArray<StorageIndex> constexpr invtrans(StorageIndex const& aQfrom,
                                                 ScalarType const& aEvent) const
      noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.invtrans_impl(aQfrom, aEvent);
    }

    /*! \brief Returns EventsSet relative to state q
     *
     * @param aQ A state on the sys
     */
    EventsSet<NEvents> constexpr getStateEvents(StorageIndex const& aQ) const
      noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.getStateEvents_impl(aQ);
    }

    /*! \brief Returns EventsSet relative to state inv q
     *
     * @param aQ A state on the sys
     */
    EventsSet<NEvents> constexpr getInvStateEvents(StorageIndex const& aQ) const
      noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.getInvStateEvents_impl(aQ);
    }

    /*! \brief Invert graph
     *
     */
    void constexpr allocateInvertedGraph() const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.allocateInvertedGraph_impl();
    }

    /*! \brief Free inverted graph
     *
     */
    void constexpr clearInvertedGraph() const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.clearInvertedGraph_impl();
    }

    /*! \brief Observer property checker
     *
     * @param aLang Language
     */
    bool constexpr checkObsProp(EventsSet<NEvents> const& aAlphabet) const
      noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys.checkObsProp_impl(aAlphabet);
    }

    /*! \brief Search observer property
     *
     * @param aLang Language
     */
    bool constexpr searchObsProp(EventsSet<NEvents> const&) const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys;
    }

    /*! \brief Projection operation
     *
     * @param aAlphanet Events set which the system will be projected.
     */
    RealDESystem constexpr& proj(EventsSet<NEvents> const& aAlphabet) noexcept
    {
        RealDESystem& sys = static_cast<RealDESystem&>(*this);
        return sys.proj_impl(aAlphabet);
    }

    /*! \brief Inverse projection operation
     *
     * @param aLang Set of events
     */
    RealDESystem constexpr& invproj(EventsSet<NEvents> const&) const noexcept
    {
        RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
        return sys;
    }

protected:
    /*! \brief Current system's states number
     *
     *
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

private:
    /*! \brief Derived class is a friend
     *
     * It avoids some issues with CRTP pattern. Constructors should be private,
     * so only friend classes, aka derived classes, are allowed to call them.
     */
    friend RealDESystem;

    /*! \brief DESystem constructor
     * \details Create a system with the defined params.
     *
     * @param aStatesNumber Number of states of the system
     * @param aInitState System's initial state
     */
    DESystemBase(StorageIndex const& aStatesNumber,
                 StorageIndex const& aInitState);

    /*! \brief Default constructor
     * \details Creates an empty system with 0 states and initial
     * states 0. Marked states is undefined.
     *
     */
    DESystemBase();

    /*! \brief DESystem destructor
     * \details It is virtual for avoiding issues with polymorphism and memory
     * leaks.
     */
    virtual ~DESystemBase() = default;

    /*! \brief Move constructor
     * \details Enable move semantics.
     */
    DESystemBase(DESystemBase&&) = default;
};
}

// include method definitions
#include "cldes/src/des/DESystemBaseCore.hpp"

#endif // DESYSTEMBASE_HPP
