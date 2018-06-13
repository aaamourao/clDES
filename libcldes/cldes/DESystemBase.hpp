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

/*
 * Forward declarion of DESystemBase's friends class TransitionProxy. A
 * transition is an element of the adjascency matrix which implements
 * the des graph.
 */
template<uint8_t NEvents, typename StorageIndex>
class TransitionProxy;

/*! \brief Discrete-Events System Base
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type use for indexing the ajacency matrix
 */
template<uint8_t NEvents, typename StorageIndex>
class DESystemBase
{
public:
    /*! \brief Adjacency matrix of bitarrays implementing a graph
     *
     * The graph represents the DES automata:
     * Non zero element: contains at least one transition
     * Each non 0 bit of each element: event that lead to the next stage
     *
     * row index: from state
     * col index: to state
     */
    using GraphHostData =
      Eigen::SparseMatrix<EventsSet<NEvents>, Eigen::RowMajor>;

    /*! \brief Set of states type
     */
    using StatesSet = std::set<StorageIndex>;

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
                                 StorageIndex const& aInitState,
                                 StatesSet& aMarkedStates)
    {
        states_number_ = aStatesNumber;
        init_state_ = aInitState;
        marked_states_ = aMarkedStates;

        // Resize graphs and do not preserve elements
        graph_.resize(states_number_, states_number_);

        // Change graphs storage type to CSR
        graph_.makeCompressed();

        // Reserve mem for hash tables for avoiding re-hashing
        states_events_ = StatesEventsTable(states_number_);
        inv_states_events_ = StatesEventsTable(states_number_);
    };

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

    /*! \brief DESystem destructor
     */
    virtual ~DESystemBase() = default;

    /*! \brief Graph getter
     */
    inline GraphHostData GetGraph() const { return graph_; };

    /*! \brief Returns state set containing the accessible part of automa
     */
    virtual StatesSet AccessiblePart() = 0;

    /*! \brief Returns state set containing the coaccessible part of automata
     */
    virtual StatesSet CoaccessiblePart() = 0;

    /*! \brief Returns States Set which is the Trim part of the system
     */
    virtual StatesSet TrimStates() = 0;

    /*! \brief Returns DES which is the Trim part of this
     */
    virtual void Trim() = 0;

    /*! \brief Returns value of the specified transition
     *
     * @param aLin Element's line
     * @param aCol Element's column
     */
    inline EventsSet<NEvents> const operator()(StorageIndex const& aLin,
                                               StorageIndex const& aCol) const
    {
        return graph_.coeff(aLin, aCol);
    }

    /*! \brief Returns value of the specified transition
     *
     * @param aLin Element's line
     * @param aCol Element's column
     */
    inline TransitionProxy<NEvents, StorageIndex> operator()(
      StorageIndex const& aLin,
      StorageIndex const& aCol)
    {
        return TransitionProxy<NEvents, StorageIndex>(this, aLin, aCol);
    }

    /*! \brief Returns number of states of the system
     *
     * Returns states_value_ by value.
     */
    StorageIndex Size() const { return states_number_; }

    /*! \brief Insert events
     *
     * @params aEvents Set containing all the new events of the current system
     */
    virtual void InsertEvents(EventsSet<NEvents> const& aEvents) = 0;

    /*
     * TODO:
     * getters
     * setters: e.g. remove transition, remove state
     * ...
     */
protected:
    // Proxy to a matrix element
    friend class TransitionProxy<NEvents, StorageIndex>;

    /*! \brief Default constructor disabled
     *
     * Declare default constructor as protected to avoid the class user of
     * calling it.
     */
    DESystemBase() = default;

    /*! \brief Let the derived class know that a new element was inserted
     *
     * @param aQfrom State from
     * @param aQto State to
     */
    virtual void NewTransition_(StorageIndex const& aQfrom,
                                StorageIndex const& aQto) = 0;

    /*! \brief Graph represented by an adjascency matrix
     *
     * A sparse matrix which represents the automata as a graph in an
     * adjascency matrix. It is implemented as a CSR scheme.
     *
     * Non zero element: transition from <row index> to <col index> when events
     * represented by the set bit indexes occurs.
     *
     * e.g. M(2, 3) = 101; Transition from state 2 to state 3 with the condition
     * event 0 OR event 2.
     */
    GraphHostData graph_;

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

    /*! \brief Current system's marked states
     *
     * Hold all marked states. Cannot be const, since the automata can be
     * cut, and some marked states may be deleted.
     */
    StatesSet marked_states_;

    /*! \brief System's events hash table
     *
     * A hash table containing all the events that matter for the current
     * system.
     */
    EventsSet<NEvents> events_;

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

/*! \brief Alias for graph 3-tuple
 *
 * (s_from, s_to, transition_events)
 */
template<uint8_t NEvents>
using Triplet = Eigen::Triplet<EventsSet<NEvents>>;
}

// include transition proxy class
#include "cldes/TransitionProxy.hpp"

#endif // DESYSTEMBASE_HPP
