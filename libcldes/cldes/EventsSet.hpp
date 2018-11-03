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
#include <istream>
#include <memory>
#include <ostream>
#include <string>

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
    explicit constexpr EventsSet<NEvents>() noexcept
      : std::bitset<NEvents>{} {};

    explicit constexpr EventsSet<NEvents>(unsigned long const val) noexcept
      : std::bitset<NEvents>{ val }
    {}

    template<class _CharT, class _Traits, class _Alloc>
    explicit constexpr EventsSet<NEvents>(
      std::basic_string<_CharT, _Traits, _Alloc> const& __s,
      size_t __position = 0) noexcept
      : std::bitset<NEvents>{ __s, __position } {};

    template<class _CharT, class _Traits, class _Alloc>
    EventsSet<NEvents>(const std::basic_string<_CharT, _Traits, _Alloc>& __s,
                       size_t __position,
                       size_t __n)
      : std::bitset<NEvents>{ __s, __position, __n }
    {}

    template<typename _CharT>
    explicit constexpr EventsSet<NEvents>(
      _CharT const* __str,
      typename std::basic_string<_CharT>::size_type __n =
        std::basic_string<_CharT>::npos,
      _CharT __zero = _CharT('0'),
      _CharT __one = _CharT('1')) noexcept
      : std::bitset<NEvents>{ __str, __n, __zero, __one }
    {}

    constexpr EventsSet<NEvents>(std::bitset<NEvents> const& val) noexcept
      : std::bitset<NEvents>{ val }
    {}

    constexpr EventsSet<NEvents>(std::bitset<NEvents>& val) noexcept
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

    EventsSet<NEvents>& operator+=(EventsSet<NEvents> const& value) noexcept
    {
        std::bitset<NEvents>::operator|=(value);
        return *this;
    }

    constexpr operator bool() const { return this->any(); }

    bool constexpr operator!=(EventsSet<NEvents> const& value) const noexcept
    {
        return std::bitset<NEvents>::operator!=(value);
    }

    bool constexpr operator==(EventsSet<NEvents> const& value) const noexcept
    {
        return std::bitset<NEvents>::operator==(value);
    }
};

/*! \brief Overload operator+ from base class
 */
template<uint8_t NEvents>
inline EventsSet<NEvents> constexpr
operator+(EventsSet<NEvents> const& aLhs, EventsSet<NEvents> const& aRhs)
{
    return aLhs | aRhs;
}

template<uint8_t NEvents>
inline EventsSet<NEvents> const&
conj(EventsSet<NEvents> const& x) noexcept
{
    return x;
}
template<uint8_t NEvents>
inline EventsSet<NEvents> const&
real(EventsSet<NEvents> const& x)
{
    return x;
}
template<uint8_t NEvents>
inline EventsSet<NEvents> constexpr imag(EventsSet<NEvents> const&) noexcept
{
    return 0;
}
template<uint8_t NEvents>
inline EventsSet<NEvents> constexpr abs(EventsSet<NEvents> const& x) noexcept
{
    return x;
}

// It does not mean anything for an events set, but eigen requires it
// WARNING: It could lead to errors using different lin algebra algorithms
template<uint8_t NEvents>
inline EventsSet<NEvents> constexpr abs2(EventsSet<NEvents> const& x) noexcept
{
    return x;
}

// WARNING: It could lead to errors using different lin algebra algorithms
template<uint8_t NEvents>
inline EventsSet<NEvents> constexpr sqrt(EventsSet<NEvents> const& x) noexcept
{
    return x;
}
} // namespace cldes

namespace std {

template<class _CharT, class _Traits, uint8_t _Nb>
std::basic_istream<_CharT, _Traits>&
operator>>(std::basic_istream<_CharT, _Traits>& __is,
           cldes::EventsSet<_Nb>& __x)
{
    std::shared_ptr<std::bitset<_Nb>> baseslice =
      std::dynamic_pointer_cast<std::bitset<_Nb>>(__x);
    return std::operator>>(__is, *baseslice);
}

template<class _CharT, class _Traits, uint8_t _Nb>
std::basic_ostream<_CharT, _Traits>&
operator<<(std::basic_ostream<_CharT, _Traits>& __os,
           cldes::EventsSet<_Nb> const& __x)
{
    std::shared_ptr<std::bitset<_Nb>> baseslice =
      std::dynamic_pointer_cast<std::bitset<_Nb>>(__x);
    return std::operator<<(__os, *baseslice);
}

} // namespace std

// Add eigen support to EventsSet scalar type
namespace Eigen {
template<uint8_t NEvents>
struct NumTraits<cldes::EventsSet<NEvents>>
  : GenericNumTraits<cldes::EventsSet<NEvents>>
{
    typedef cldes::EventsSet<NEvents> Real;
    typedef cldes::EventsSet<NEvents> NonInteger;
    typedef cldes::EventsSet<NEvents> Nested;

    enum
    {
        IsComplex = 0,
        IsInteger = 1,
        IsSigned = 0,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 1,
        MulCost = 1
    };
};
}

#endif // CLDES_EVENTSSET_HPP
