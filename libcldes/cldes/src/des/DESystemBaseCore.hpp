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
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::setEvents(
  EventsSet<NEvents> const& aEvents) noexcept
{
    events_ = aEvents;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::setStatesNumber(
  StorageIndex const& aStNum) noexcept
{
    states_number_ = aStNum;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::setInitialState(
  StorageIndex const& aInitState) noexcept
{
    init_state_ = aInitState;
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::insertMarkedState(
  StorageIndex const& aSt) noexcept
{
    marked_states_.emplace(aSt);
}

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
void
DESystemBase<NEvents, StorageIndex, RealDESystem>::setMarkedStates(
  StatesSet const& aStSet) noexcept
{
    marked_states_ = aStSet;
}

}
