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

 File: cldes/Constants.hpp
 Description: Library constants and default alias definitions
 =========================================================================
*/
#ifndef CLDES_CONSTANTS_HPP
#define CLDES_CONSTANTS_HPP

// Need to define it, since it could be compiled on platforms which support
// newer OpenCL standards.
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <Eigen/Sparse>
#include <bitset>
#include <limits>

namespace cldes {

// Host adjascency matrix base type which represents an array of bits
using ScalarType = uint8_t;

/*
 * Global const definitions: preferred over defines
 */
// Max number of events
ScalarType const kMaxEvents = std::numeric_limits<ScalarType>::max();

// Max number of events on GPU objects
// TODO: Change this limit to an OpenCL definition for different platforms
ScalarType const kMaxEventsGPU = std::numeric_limits<uint64_t>::max();

} // namespace cldes

#endif // CLDES_CONSTANTS_HPP
