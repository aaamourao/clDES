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

 File: transition_proxy.hpp
 Description: Proxy to an element of the graph_ data member from DESystem.
 =========================================================================
*/

#ifndef TRANSITION_PROXY_HPP
#define TRANSITION_PROXY_HPP

#include "constants.hpp"

namespace cldes {

class DESystem;

class TransitionProxy {
public:
    TransitionProxy(DESystem *const aSysPtr, cldes_size_t const &aLin,
                    cldes_size_t const &aCol);
    TransitionProxy &operator=(ScalarType aTransitionValue);
    operator ScalarType();

private:
    TransitionProxy();
    DESystem * sys_ptr_;
    cldes_size_t const lin_;
    cldes_size_t const col_;
};

}  // namespace cldes

#endif  // TRANSITION_PROXY_HPP
