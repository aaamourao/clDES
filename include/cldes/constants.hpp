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

 File: constants.hpp
 Description: Library constants
 =========================================================================
*/
#ifndef CLDES_CONSTANTS_HPP
#define CLDES_CONSTANTS_HPP

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

//#include <CL/cl.hpp>
#include <bitset>
#include <Eigen/Sparse>

namespace cldes {
// Host adjascency matrix base type which represents an array of bits
using ScalarType = size_t;

template<size_t NEvents, typename StorageIndex>
using Triplet = Eigen::Triplet<std::bitset<NEvents>,
                               typename std::make_signed<StorageIndex>::type>;

template<typename StorageIndex>
using BitTriplet =
  Eigen::Triplet<bool, typename std::make_signed<StorageIndex>::type>;
} // namespace cldes

namespace std {

template<std::size_t size>
inline bitset<size>
operator+(bitset<size> const& rhs, bitset<size> const& lhs)
{
    return rhs | lhs;
}

} // namespace std

#endif // CLDES_CONSTANTS_HPP
