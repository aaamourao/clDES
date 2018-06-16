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

 File: cldes/TransitionProxy.hpp
 Description: Proxy to an element of the graph_ data member from
 DESystem.
 =========================================================================
*/

#ifndef TRANSITION_PROXY_HPP
#define TRANSITION_PROXY_HPP

#include "cldes/Constants.hpp"

namespace cldes {

/*! \brief Proxy to a transition of a Discrete-Events System
 *
 * @param NEvents Number of events
 * @param StorageIndex Unsigned type for storing the indexes of each state
 */
template<uint8_t NEvents, typename StorageIndex>
class TransitionProxy
{
public:
    /*! \brief TransitionProxy unique constructor
     *
     * Feed const data-members.
     *
     * @param aSysPtr Raw pointer to owner of the elem which will be proxied
     * @param aLin Line of the element
     * @param aCol Column of the element
     */
    TransitionProxy(DESystem<NEvents, StorageIndex>* const aSysPtr,
                    StorageIndex const& aLin,
                    StorageIndex const& aCol);

    /*! \brief Override operator "=" from TransitionProxy class
     *
     * Override operator "=" for tracking when the device graph from the
     * related DESystem object is outdated.
     *
     * @param aTransitionValue graph_(lin_, col_) new value
     */
    TransitionProxy& operator=(ScalarType aTransitionValue);

    /*! \brief Override cast to ScalarType
     *
     * Operator used when element from graph_ is accessed, but its value is not
     * changed.
     */
    operator EventsSet<NEvents>() const;

protected:
    /*! \brief Disabled default constructor
     *
     * There is no use for the default constructor.
     */
    TransitionProxy();

private:
    /*! \brief Raw pointer to DESystemBase object
     *
     * Raw pointer to the owner of the proxied element.
     */
    DESystem<NEvents, StorageIndex>* sys_ptr_;

    /*! \brief Element line.
     *
     * Line where the elem is located.
     */
    StorageIndex const lin_;

    /*! \brief Column line.
     *
     * Column where the elem is located.
     */
    StorageIndex const col_;
};

} // namespace cldes

// including implementation
#include "cldes/src/des/TransitionsProxyCore.hpp"

#endif // TRANSITION_PROXY_HPP
