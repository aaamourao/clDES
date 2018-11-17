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

    StorageIndex const sparcitypattern = aproxSpacPat_(aVirtualSys);
    aVirtualSys.resizeStatesEvents(aVirtualSys.states_number_);
    aVirtualSys.trans_number_ = 0;
    aVirtualSys.triplet_.reserve(sparcitypattern);
    for (StorageIndex qfrom = 0; qfrom < aVirtualSys.states_number_; ++qfrom) {
        aVirtualSys.setStateEvents(qfrom, aVirtualSys.getStateEvents(qfrom));
        aVirtualSys.setInvStateEvents(qfrom,
                                      aVirtualSys.getInvStateEvents(qfrom));
        auto event = 0u;
        auto event_it = aVirtualSys.getStateEvents(qfrom);
        while (event_it != 0) {
            if (event_it.test(0)) {
                auto const qto = aVirtualSys.trans(qfrom, event);
                aVirtualSys.triplet_.push_back(Triplet<NEvents>(
                  qfrom, qto, EventsSet<NEvents>{ 1ul << event }));
                ++aVirtualSys.trans_number_;
            }
            ++event;
            event_it >>= 1;
        }
    }
    return;
}

template<class SysT_l, class SysT_r>
void
synchronizeStage2(SyncSysProxy<SysT_l, SysT_r>& aVirtualSys) noexcept
{
    using StorageIndex = typename SysTraits<SysT_l>::Si_;

    StorageIndex const sparcitypattern = aproxSpacPat_(aVirtualSys);
    SparseStatesMap<StorageIndex> statesmap;
    aVirtualSys.setStatesNumber(aVirtualSys.virtual_states_.size());

    // TODO: It SHOULD be returned in the future in a pair
    StorageIndex cst = 0;
    for (StorageIndex s : aVirtualSys.virtual_states_) {
        statesmap[s] = cst;
        ++cst;
    }
    aVirtualSys.virtual_states_.clear();
    for (StorageIndex s0 : aVirtualSys.sys0_.getMarkedStates()) {
        for (StorageIndex s1 : aVirtualSys.sys1_.getMarkedStates()) {
            StorageIndex const key = s1 * aVirtualSys.n_states_sys0_ + s0;
            if (statesmap.contains(key)) {
                aVirtualSys.insertMarkedState(statesmap[key]);
            }
        }
    }
    processVirtSys_(aVirtualSys, sparcitypattern, std::move(statesmap));
    return;
}

template<class SysT_l, class SysT_r>
#ifdef __clang__
inline long unsigned
#elif __GNUC__
long unsigned
#endif
aproxSpacPat_(SyncSysProxy<SysT_l, SysT_r> const& aV) noexcept
{
    if (aV.trans_number_ > 0) {
        return aV.trans_number_;
    }
    return (aV.events_.count() * aV.states_number_) / 3;
}

template<class SysT_l, class SysT_r>
#ifdef __clang__
inline void
#elif __GNUC__
void
#endif
processVirtSys_(SyncSysProxy<SysT_l, SysT_r>& aVirtualSys,
                unsigned long const& aSparcityPattern,
                SparseStatesMap_t<SysT_l>&& aStatesMap) noexcept
{
    uint8_t static constexpr NEvents = SysTraits<SysT_l>::Ne_;

    aVirtualSys.triplet_.reserve(aSparcityPattern);
    aVirtualSys.trans_number_ = 0;
    while (!aVirtualSys.transtriplet_.empty()) {
        auto q_trans = aVirtualSys.transtriplet_.back();
        auto const q = q_trans.first;
        while (!q_trans.second->empty()) {
            auto const qto_e = q_trans.second->back();
            auto const qto = qto_e.first;
            if (aStatesMap.contains(qto)) {
                aVirtualSys.triplet_.push_back(
                  Triplet<NEvents>(aStatesMap[q],
                                   aStatesMap[qto],
                                   EventsSet<NEvents>{ 1ul << qto_e.second }));
                ++aVirtualSys.trans_number_;
            }
            q_trans.second->pop_back();
        }
        delete q_trans.second;
        aVirtualSys.transtriplet_.pop_back();
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

template<class SysT_l, class SysT_r, typename StorageIndex>
inline void
removeBadStates_(SyncSysProxy<SysT_l, SysT_r> const& aVirtualSys,
                 transMap_t<SysT_l>& aC,
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
                            delete aC[s];
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
    uint8_t static constexpr NEvents = SysTraits<SysT_l>::Ne_;
    using StorageIndex = typename SysTraits<SysT_l>::Si_;

    SyncSysProxy<SysT_l, SysT_r> virtualsys{ aP, aE };
    EventsSet<NEvents> non_contr_bit;
    EventsSet<NEvents> p_non_contr_bit;

    for (cldes::ScalarType event : aNonContr) {
        if (aP.getEvents().test(event)) {
            p_non_contr_bit.set(event);
            if (virtualsys.getEvents().test(event)) {
                non_contr_bit.set(event);
            }
        }
    }
    transMap<StorageIndex>&& c = computeSupCStates_(
      virtualsys, std::move(non_contr_bit), std::move(p_non_contr_bit), aP);
    virtualsys.virtual_states_.reserve(c.size());
    virtualsys.transtriplet_.reserve(c.size());
    for (auto tr : c) {
        virtualsys.virtual_states_.push_back(tr.first);
        virtualsys.transtriplet_.push_back(tr);
        virtualsys.trans_number_ += (*tr.second).size() / 2;
    }
    c.clear();

    auto sys = std::move(DESystem<NEvents, StorageIndex>(virtualsys));
    sys.trim();

    return sys;
}

template<class SysT_l, class SysT_r>
#ifdef __clang__
inline transMap_t<SysT_l>
#elif __GNUC__
transMap_t<SysT_l>
#endif
computeSupCStates_(SyncSysProxy<SysT_l, SysT_r> const& aVirtualSys,
                   EventsSet_t<SysT_l> const&& aNonContrBit,
                   EventsSet_t<SysT_l> const&& aPNonContrBit,
                   SysT_l const& aP) noexcept
{
    using StorageIndex = typename SysTraits<SysT_l>::Si_;
    transMap<StorageIndex> c;
    StatesTableHost<StorageIndex> rmtable;
    StatesStack<StorageIndex> f;
    f.push(aVirtualSys.getInitialState());
    aVirtualSys.allocateInvertedGraph();
    while (!f.empty()) {
        auto const q = f.top();
        f.pop();
        if (!rmtable.contains(q) && !c.contains(q)) {
            auto const q_events = aVirtualSys.getStateEvents(q);
            auto const in_ncqx =
              aPNonContrBit & aP.getStateEvents(q % aVirtualSys.n_states_sys0_);
            if ((in_ncqx & q_events) != in_ncqx) {
                removeBadStates_(aVirtualSys, c, q, aNonContrBit, rmtable);
            } else {
                c[q] = new InvArgtrans<StorageIndex>();
                cldes::ScalarType event = 0;
                auto event_it = q_events;
                while (event_it.any()) {
                    if (event_it.test(0)) {
                        auto const fsqe = aVirtualSys.trans(q, event);
                        if (!rmtable.contains(fsqe)) {
                            if (!c.contains(fsqe)) {
                                f.push(fsqe);
                            }
                        }
                        c[q]->push_back(std::make_pair(fsqe, event));
                    }
                    ++event;
                    event_it >>= 1;
                }
            }
        }
    }
    aVirtualSys.clearInvertedGraph();
    return c;
}

// template<class SysT_l, class SysT_r>
// GenericSystem
// createBinExprTree(DESVector<NEvents, StorageIndex> const& aSystems)
// {
//     std::vector<std::shared_ptr<DESystemBase>> sys;
//     std::vector<std::shared_ptr<DESystemBase>> nodes_ref;
//
//     for (auto s : aSystems) { // initialize tree
//         auto node = std::make_shared<DESystem>(s);
//         sys.push_back(node);
//         nodes_ref.push_back(node);
//     }
//     while (sys.size() != 1) {
//         auto cp_sys = std::move(sys);
//         if (cp_sys.size() % 2 != 0) {
//             std::shared_ptr<DESystemBase> node =
//               cp_sys.back(); // node is a shared ptr
//             cp_sys.pop_back();
//             sys.push_back(node);
//         }
//         size_t processed_nodes = 0ul;
//         while (!cp_sys.empty()) { // So it has even number of items
//             std::shared_ptr<DESystemBase> lhs = cp_sys.back();
//             cp_sys.pop_back();
//             std::shared_ptr<DESystemBase> rhs = cp_sys.back();
//             cp_sys.pop_back();
//             std::shared_ptr<DESystemBase> node =
//               std::make_shared<SyncSysProxy>(SyncSysProxy{ *lhs, *rhs });
//             nodes_ref.push_back(node);
//             sys.push_back(node);
//             processed_nodes += 2;
//         }
//     }
//     return std::make_pair(sys[0], nodes_ref);
// }

// template<uint8_t NEvents, typename StorageIndex>
// DESystem<NEvents, StorageIndex>
// supC(DESVector<NEvents, StorageIndex> const& aPlants,
//                 DESVector<NEvents, StorageIndex> const& aSpecs,
//                 EventsTableHost const& aNonContr)
// {
//     using BinExprTree = BinExprTree<NEvents, StorageIndex>;
//
//     BinExprTree plant = GenBinExprTree(aPlants);
//     BinExprTree spec = GenBinExprTree(aSpecs);
//
//     auto const supervisor = supC<NEvents, StorageIndex>(
//       *(plant.first), *(spec.first), aNonContr);
//
//     return supervisor;
// }

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
