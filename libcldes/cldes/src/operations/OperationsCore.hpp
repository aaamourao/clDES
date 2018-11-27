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

 File: cldes/operations/Operations.hpp
 Description: Definition of operation functions: Parallel composition,
 virtual parallel composition and supervisor synthesis.
 =========================================================================
*/

#ifndef OPERATIONS_CORE_HPP
#define OPERATIONS_CORE_HPP

namespace cldes {
namespace op {

// TODO: put all stage2 on impl and make a single function
template<class SysT_l, class SysT_r>
void
synchronizeEmptyStage2(SyncSysProxy<SysT_l, SysT_r>& aVirtualSys) noexcept
{
    uint8_t static constexpr NEvents = SysTraits<SysT_l>::Ne_;
    using StorageIndex = typename SysTraits<SysT_l>::Si_;

    // aprox
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_ / 3;
    aVirtualSys.resizeStatesEvents(aVirtualSys.states_number_);
    aVirtualSys.trans_number_ = 0;
    aVirtualSys.triplet_.reserve(sparcitypattern);
#ifdef CLDES_OPENMP_ENABLED
#pragma omp parallel
    {
        std::vector<Triplet<NEvents>> triplet_parallel;
#pragma omp for nowait
        for (StorageIndex qfrom = 0; qfrom < aVirtualSys.states_number_;
             ++qfrom) {
#else
    {
        for (StorageIndex qfrom = 0; qfrom < aVirtualSys.states_number_;
             ++qfrom) {
#endif
            aVirtualSys.setStateEvents(qfrom,
                                       aVirtualSys.getStateEvents(qfrom));
            aVirtualSys.setInvStateEvents(qfrom,
                                          aVirtualSys.getInvStateEvents(qfrom));
            auto event = 0u;
            auto event_it = aVirtualSys.getStateEvents(qfrom);
            while (event_it != 0) {
                if (event_it.test(0)) {
                    auto const qto = aVirtualSys.trans(qfrom, event);
#ifdef CLDES_OPENMP_ENABLED
                    triplet_parallel.push_back(
#else
                    aVirtualSys.triplet_.push_back(
#endif
                      Triplet<NEvents>(
                        qfrom, qto, EventsSet<NEvents>{ 1ul << event }));
                    ++aVirtualSys.trans_number_;
                }
                ++event;
                event_it >>= 1;
            }
        }
#ifdef CLDES_OPENMP_ENABLED
#pragma omp critical
        aVirtualSys.triplet_.insert(aVirtualSys.triplet_.end(),
                                    triplet_parallel.begin(),
                                    triplet_parallel.end());
#endif
    }
    return;
}

template<class SysT_l, class SysT_r, class StTabT, typename StorageIndex>
inline void
removeBadStates_(SyncSysProxy<SysT_l, SysT_r> const& aVirtualSys,
                 StTabT& aC,
                 StorageIndex const& aQ,
                 EventsSet_t<SysT_l> const& aNonContrBit,
                 StatesTableHost_t<SysT_l>& aRmTable) noexcept
{
    StatesStack<StorageIndex> f;
    f.push(aQ);
    aRmTable.insert(aQ);
    while (!f.empty()) {
        auto const x = f.top();
        f.pop();
        auto q_events = aVirtualSys.getInvStateEvents(x);
        q_events &= aNonContrBit;
        cldes::ScalarType event = 0;
        while (q_events.any()) {
            if (q_events.test(0)) {
                StatesArray<StorageIndex> const finv =
                  aVirtualSys.invtrans(x, event);
                for (StorageIndex s : finv) {
                    if (!aRmTable.contains(s)) {
                        f.push(s);
                        aRmTable.insert(s);
                        if (aC.contains(s)) {
                            aC.erase(s);
                        }
                    }
                }
            }
            ++event;
            q_events >>= 1;
        }
    }
    return;
}

template<class SysT_l, class SysT_r>
DESystem_t<SysT_l>
supC(SysT_l const& aP,
     SysT_r const& aE,
     EventsTableHost const& aNonContr) noexcept
{
    uint8_t constexpr NEvents = SysTraits<SysT_l>::Ne_;
    using StorageIndex = typename SysTraits<SysT_l>::Si_;

    DESystem<NEvents, StorageIndex> sys = DESystem<NEvents, StorageIndex>(
      SuperProxy<SysT_l, SysT_r>{ aP, aE, aNonContr });

    return sys;
}

template<class SysT>
SysT&
proj(SysT const& aSys, EventsSet_t<SysT> const&) noexcept
{
    return aSys;
}

template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
proj(DESystem<NEvents, StorageIndex> const& aSys,
     EventsSet<NEvents> const&) noexcept
{
    return aSys;
}

} // namespace op
} // namespace cldes

#endif // OPERATIONS_CORE_HPP
