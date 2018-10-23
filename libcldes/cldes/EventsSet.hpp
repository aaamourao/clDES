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

#include <Eigen/Core>
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
class EventsSet : public std::bitset<NEvents>
{
public:
    explicit EventsSet<NEvents>()
      : std::bitset<NEvents>{} {};

    explicit EventsSet<NEvents>(unsigned long const val)
      : std::bitset<NEvents>{ val }
    {}

    EventsSet<NEvents>(std::bitset<NEvents> const& val)
      : std::bitset<NEvents>{ val }
    {}
    EventsSet<NEvents>(std::bitset<NEvents>& val)
      : std::bitset<NEvents>{ val }
    {}

    EventsSet<NEvents>(EventsSet<NEvents>&) = default;
    EventsSet<NEvents>(EventsSet<NEvents> const&) = default;
    EventsSet<NEvents>(EventsSet<NEvents>&&) = default;
    ~EventsSet<NEvents>() = default;

    EventsSet<NEvents>& operator=(EventsSet<NEvents>&) = default;
    EventsSet<NEvents>& operator=(unsigned long const val)
    {
        std::bitset<NEvents>::operator=(val);
        return *this;
    }
    EventsSet<NEvents>& operator=(EventsSet<NEvents> const&) = default;
    EventsSet<NEvents>& operator=(EventsSet<NEvents>&&) = default;

    EventsSet<NEvents>& operator+=(EventsSet<NEvents> const& value)
    {
        std::bitset<NEvents>::operator|=(value);
        return *this;
    }

    operator bool() const { return this->any(); }

    bool operator!=(EventsSet<NEvents> const& value)
    {
        return std::bitset<NEvents>::operator!=(value);
    }

    bool operator==(EventsSet<NEvents> const& value)
    {
        return std::bitset<NEvents>::operator==(value);
    }
};

/*! \brief Overload operator+ from base class
 */
template<uint8_t NEvents>
inline EventsSet<NEvents>
operator+(EventsSet<NEvents> const& aLhs, EventsSet<NEvents> const& aRhs)
{
    return aLhs | aRhs;
}

template<uint8_t NEvents>
inline const EventsSet<NEvents>&
conj(const EventsSet<NEvents>& x)
{
    return x;
}
template<uint8_t NEvents>
inline const EventsSet<NEvents>&
real(const EventsSet<NEvents>& x)
{
    return x;
}
template<uint8_t NEvents>
inline EventsSet<NEvents>
imag(const EventsSet<NEvents>&)
{
    return 0;
}
template<uint8_t NEvents>
inline EventsSet<NEvents>
abs(const EventsSet<NEvents>& x)
{
    return x;
}
template<uint8_t NEvents>
inline EventsSet<NEvents>
abs2(const EventsSet<NEvents>& x)
{
    return x;
}

template<uint8_t NEvents>
inline EventsSet<NEvents>
sqrt(const EventsSet<NEvents>& x)
{
    return x;
}
} // namespace cldes

// Add eigen support to EventsSet scalar type
namespace Eigen {
template<uint8_t NEvents>
struct NumTraits<cldes::EventsSet<NEvents>>
  : GenericNumTraits<cldes::EventsSet<NEvents>>
{
    typedef cldes::EventsSet<NEvents> Real;
    typedef cldes::EventsSet<NEvents> Integer;
    typedef cldes::EventsSet<NEvents> Nested;

    enum
    {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 0,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 1,
        MulCost = 1
    };
};
}

#endif // CLDES_EVENTSSET_HPP
