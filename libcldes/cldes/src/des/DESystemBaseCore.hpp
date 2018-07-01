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

 File: cldes/src/des/DESystemBaseCore.hpp
 Description: DESystemBase abstract class methods definition (the
 non-abstract methods, of course)
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
namespace cldes {
template<uint8_t NEvents, typename StorageIndex>
DESystemBase<NEvents, StorageIndex>::DESystemBase(
  StorageIndex const& aStatesNumber,
  StorageIndex const& aInitState)
{
    states_number_ = aStatesNumber;
    init_state_ = aInitState;
}

template<uint8_t NEvents, typename StorageIndex>
DESystemBase<NEvents, StorageIndex>::DESystemBase()
{
    states_number_ = 0;
    init_state_ = 0;
    events_ = EventsSet<NEvents>{};
    marked_states_ = StatesSet{};
}

template<uint8_t NEvents, typename StorageIndex>
StorageIndex
DESystemBase<NEvents, StorageIndex>::Size() const
{
    return states_number_;
}

template<uint8_t NEvents, typename StorageIndex>
EventsSet<NEvents>
DESystemBase<NEvents, StorageIndex>::GetEvents() const
{
    return events_;
}

template<uint8_t NEvents, typename StorageIndex>
StorageIndex
DESystemBase<NEvents, StorageIndex>::GetStatesNumber() const
{
    return states_number_;
}

template<uint8_t NEvents, typename StorageIndex>
StorageIndex
DESystemBase<NEvents, StorageIndex>::GetInitialState() const
{
    return init_state_;
}

template<uint8_t NEvents, typename StorageIndex>
typename DESystemBase<NEvents, StorageIndex>::StatesSet
DESystemBase<NEvents, StorageIndex>::GetMarkedStates() const
{
    return marked_states_;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetEvents(
  EventsSet<NEvents> const& aEvents)
{
    events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetStatesNumber(StorageIndex const& aStNum)
{
    states_number_ = aStNum;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetInitialState(
  StorageIndex const& aInitState)
{
    init_state_ = aInitState;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::InsertMarkedState(StorageIndex const& aSt)
{
    marked_states_.emplace(aSt);
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetMarkedStates(StatesSet const& aStSet)
{
    marked_states_ = aStSet;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::ResizeStatesEvents(
  StorageIndex const& aSize)
{
    states_events_.resize(aSize);
    inv_states_events_.resize(aSize);
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetStatesEvents(
  StatesEventsTable const& aEvents)
{
    states_events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetInvStatesEvents(
  StatesEventsTable const& aEvents)
{
    inv_states_events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetStateEvents(
  StorageIndex const& aQ,
  EventsSet<NEvents> const& aEvent)
{
    states_events_[aQ] = aEvent;
}

template<uint8_t NEvents, typename StorageIndex>
void
DESystemBase<NEvents, StorageIndex>::SetInvStateEvents(
  StorageIndex const& aQ,
  EventsSet<NEvents> const& aEvent)
{
    inv_states_events_[aQ] = aEvent;
}
}
