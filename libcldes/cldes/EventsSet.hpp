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

 File: cldes/EventsSet.hpp
 Description: Events set encapsulated as a std::bitset with overloaded
 operators.
 =========================================================================
*/
#ifndef CLDES_EVENTSSET_HPP
#define CLDES_EVENTSSET_HPP

#include <bitset>

namespace cldes {

/*! \brief Bit array representing an events set
 *
 * So far, it could be alias. It is more efficient than defining a class for
 * overloading a single operator.
 *
 * However, it leads to an issue when deducing template parameters on
 * cldes::op::RemoveBadStates and static cldes::op::__TransitionVirtualInv
 *
 * Each bit represent a different event.
 * 0 -> does not contain event
 * 1 -> contains event
 * index -> Event value
 *
 * @param NEvents Number of events: max = 255
 */
template<uint8_t NEvents>
using EventsSet = std::bitset<NEvents>;

} // namespace cldes

namespace std {

/*! \brief Overload operator+ from base class
 */
template<std::size_t size>
inline bitset<size>
operator+(bitset<size> const& aLhs, bitset<size> const& aRhs)
{
    return aLhs | aRhs;
}
}

#endif // CLDES_EVENTSSET_HPP
