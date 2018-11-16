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

 File: cldes/src/operations/SuperProxyCore.hpp
 Description: includes and alias.
 =========================================================================
*/

#include <algorithm>

namespace cldes {
namespace op {

template<class SysT_l, class SysT_r>
void
synchronizeStage2(SyncSysProxy<SysT_l, SysT_r>& aVirtualSys) noexcept;

template<class SysT_l, class SysT_r>
void
synchronizeEmptyStage2(SyncSysProxy<SysT_l, SysT_r>& aVirtualSys) noexcept;

using EventsTableHost = spp::sparse_hash_set<uint8_t>;

template<class SysT_l, class SysT_r>
DESystem<SysTraits<SysT_l>::Ne_, typename SysTraits<SysT_l>::Si_>
supC(SysT_l const& aP,
     SysT_r const& aE,
     EventsTableHost const& aNonContr) noexcept;
}
}
