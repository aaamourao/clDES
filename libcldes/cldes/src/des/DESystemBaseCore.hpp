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
template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
DESystemBase<NEvents, StorageIndex, RealDESystem>::DESystemBase(
  StorageIndex const& aStatesNumber,
  StorageIndex const& aInitState)
{
    states_number_ = aStatesNumber;
    init_state_ = aInitState;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
DESystemBase<NEvents, StorageIndex, RealDESystem>::DESystemBase()
{
    states_number_ = 0;
    init_state_ = 0;
    events_ = EventsSet<NEvents>{};
    marked_states_ = StatesSet{};
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
StorageIndex
DESystemBase<NEvents, StorageIndex, RealDESystem>::Size() const
{
    return states_number_;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
EventsSet<NEvents>
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetEvents() const
{
    return events_;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
StorageIndex
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetStatesNumber() const
{
    return states_number_;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
StorageIndex
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetInitialState() const
{
    return init_state_;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
typename DESystemBase<NEvents, StorageIndex, RealDESystem>::StatesSet
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetMarkedStates() const
{
    return marked_states_;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetEvents(
  EventsSet<NEvents> const& aEvents)
{
    events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetStatesNumber(
  StorageIndex const& aStNum)
{
    states_number_ = aStNum;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetInitialState(
  StorageIndex const& aInitState)
{
    init_state_ = aInitState;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::InsertMarkedState(
  StorageIndex const& aSt)
{
    marked_states_.emplace(aSt);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetMarkedStates(
  StatesSet const& aStSet)
{
    marked_states_ = aStSet;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::ResizeStatesEvents(
  StorageIndex const& aSize)
{
    states_events_.resize(aSize);
    inv_states_events_.resize(aSize);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetStatesEvents(
  StatesEventsTable const& aEvents)
{
    states_events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetInvStatesEvents(
  StatesEventsTable const& aEvents)
{
    inv_states_events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetStateEvents(
  StorageIndex const& aQ,
  EventsSet<NEvents> const& aEvent)
{
    states_events_[aQ] = aEvent;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::SetInvStateEvents(
  StorageIndex const& aQ,
  EventsSet<NEvents> const& aEvent)
{
    inv_states_events_[aQ] = aEvent;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
bool
DESystemBase<NEvents, StorageIndex, RealDESystem>::ContainsTrans(
  StorageIndex const& aQ,
  ScalarType const& aEvents) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.containsTrans_impl(aQ, aEvents);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
bool constexpr DESystemBase<NEvents, StorageIndex, RealDESystem>::IsVirtual()
  const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.isVirtual_impl();
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
std::shared_ptr<DESystemBase<NEvents, StorageIndex, RealDESystem>>
DESystemBase<NEvents, StorageIndex, RealDESystem>::Clone() const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.clone_impl();
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
typename DESystemBase<NEvents, StorageIndex, RealDESystem>::StorageIndexSigned
DESystemBase<NEvents, StorageIndex, RealDESystem>::Trans(
  StorageIndex const& aQ,
  ScalarType const& aEvent) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.trans_impl(aQ, aEvent);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
bool
DESystemBase<NEvents, StorageIndex, RealDESystem>::ContainsInvTrans(
  StorageIndex const& aQ,
  ScalarType const& aEvent) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.containsInvTrans_impl(aQ, aEvent);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
StatesArray<StorageIndex>
DESystemBase<NEvents, StorageIndex, RealDESystem>::InvTrans(
  StorageIndex const& aQfrom,
  ScalarType const& aEvent) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.invTrans_impl(aQfrom, aEvent);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
EventsSet<NEvents>
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetStateEvents(
  StorageIndex const& aQ) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.getStateEvents_impl(aQ);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
EventsSet<NEvents>
DESystemBase<NEvents, StorageIndex, RealDESystem>::GetInvStateEvents(
  StorageIndex const& aQ) const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.getInvStateEvents_impl(aQ);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::AllocateInvertedGraph() const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.allocateInvertedGraph_impl();
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::ClearInvertedGraph() const
{
    RealDESystem const& sys = static_cast<RealDESystem const&>(*this);
    return sys.clearInvertedGraph_impl();
}
}
