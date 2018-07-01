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

 File: cldes/operations/Operations.hpp
 Description: Definition of operation functions: Parallel composition,
 virtual parallel composition and supervisor synthesis.
 =========================================================================
*/

#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "cldes/Constants.hpp"
#include "cldes/DESystem.hpp"
#include "cldes/operations/SyncSysProxy.hpp"
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <stack>
#include <tuple>

namespace cldes {

namespace op {
/*! \brief tuple representing a state of a virtual synch (stage 1)
 * \details (state_id_g0, state_id_g1)
 */
template<typename StorageIndex>
using StatesTupleHost = std::pair<StorageIndex, StorageIndex>;

/*! \brief Hash set of virtual states (stage 1)
 * \details st = state_id_g1 * g0.size() + state_id_g0
 */
template<typename StorageIndex>
using StatesTableHost = spp::sparse_hash_set<StorageIndex>;

/*! \brief Hash map type for maps a virtual state to its new index
 * \details It is necessary when states are removed.
 */
template<typename StorageIndex>
using SparseStatesMap = spp::sparse_hash_map<StorageIndex, StorageIndex>;

/*! \brief Stack of states type
 */
template<typename StorageIndex>
using StatesStack = std::stack<StorageIndex>;

/*! \brief Hash set containing events indexes
 */
using EventsTableHost = spp::sparse_hash_set<uint8_t>;

/*! \brief Calculate the parallel composition between two systems
 * \details The composed states are sorted by the right operand indexes:
 * e.g. indexes(sys0.size{3} || sys.size{2}) =
 * {0 = (0, 0), 1 = (1, 0), 2 = (2, 0), 3 = (0, 1), 4 = (1, 1), 5 = (2, 1)}
 * \warning Synchronize() is faster than the other synchronizing functions to
 * calculate the whole system. However, parallel composition between large
 * systems can occupy a lot of memory. Prefer lazy operations when this
 * is a problem
 *
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 * \return A concrete system which represents a parallel composition
 * between two systems
 */
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
Synchronize(DESystemBase<NEvents, StorageIndex> const& aSys0,
            DESystemBase<NEvents, StorageIndex> const& aSys1)
{
    DESystem<NEvents, StorageIndex> sys = DESystem<NEvents, StorageIndex>(
      SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 });

    return sys;
}

/*! \brief Lazy evaluation of the parallel composition between two systems
 * \details The composed states are sorted by the right operand indexes:
 * e.g. indexes(sys0.size{3} || sys.size{2}) =
 * {0 = (0, 0), 1 = (1, 0), 2 = (2, 0), 3 = (0, 1), 4 = (1, 1), 5 = (2, 1)}
 * \warning This funtion returns an object that evaluates the operation
 * on demand. If you need the whole system at once, use Synchronize or
 * convert with the result to DESystem using a type cast.
 *
 * @param aSys0 The left operand of the parallel composition.
 * @param aSys1 The right operand of the parallel composition.
 * \return A virtual system which represents a parallel composition
 * between two systems.
 */
template<uint8_t NEvents, typename StorageIndex>
SyncSysProxy<NEvents, StorageIndex>
SynchronizeStage1(DESystemBase<NEvents, StorageIndex> const& aSys0,
                  DESystemBase<NEvents, StorageIndex> const& aSys1)
{
    return SyncSysProxy<NEvents, StorageIndex>{ aSys0, aSys1 };
}

/*! \brief Final stage of the lazy parallel composition evaluation
 * \details Transform a virtual proxy in a concrete system. It is implicited
 * called when virtual proxy to the parallel composition operation that
 * has no value calculated is converted to a concrete system.
 *
 * TODO: Cache lazy calculated transitions when the user visit transitions and
 * use it on synchronize stage 2
 *
 * @param[out] aVirtualSys Reference to the system which will be transformed.
 * \return void
 */
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeEmptyStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    // Estimated sparcity pattern: Calculate on the same loop than statesmap
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_;

    // Reserve space for transitions
    aVirtualSys.ResizeStatesEvents(aVirtualSys.states_number_);
    aVirtualSys.triplet_.reserve(sparcitypattern);
    aVirtualSys.bittriplet_.reserve(sparcitypattern +
                                    aVirtualSys.states_number_);

    // Calculate transitions
    for (StorageIndex qfrom = 0; qfrom < aVirtualSys.states_number_; ++qfrom) {
        aVirtualSys.bittriplet_.push_back(BitTriplet(qfrom, qfrom, true));

        aVirtualSys.SetStateEvents(qfrom, aVirtualSys.GetStateEvents(qfrom));
        aVirtualSys.SetInvStateEvents(qfrom,
                                      aVirtualSys.GetInvStateEvents(qfrom));

        auto event = 0u;
        auto event_it = aVirtualSys.GetStateEvents(qfrom);
        while (event_it != 0) {
            if (event_it.test(0)) {
                auto const qto = aVirtualSys.Trans(qfrom, event);

                aVirtualSys.triplet_.push_back(Triplet<NEvents>(
                  qfrom, qto, EventsSet<NEvents>{ 1ul << event }));

                if (qfrom != static_cast<StorageIndex>(qto)) {
                    aVirtualSys.bittriplet_.push_back(
                      BitTriplet(qto, qfrom, true));
                }
            }
            ++event;
            event_it >>= 1;
        }
    }

    return;
}

/*! \brief Transform a virtual system in a real system: optmized to supervisor
 * synthesis
 *
 * Calculate transitions and other parameters of a virtual system which
 * represents a virtual parallel composition. It does not calculate the
 * states_events_ table, since supervisors are final system and may have a big
 * size
 *
 * @param[out] aVirtualSys Reference to the system which will be transformed.
 * \return void
 */
template<uint8_t NEvents, typename StorageIndex>
void
SynchronizeStage2(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys)
{
    SparseStatesMap<StorageIndex> statesmap;

    aVirtualSys.SetStatesNumber(aVirtualSys.virtual_states_.size());

    // Reserve space for transitions
    aVirtualSys.bittriplet_.reserve(aVirtualSys.states_number_);

    // Estimated sparcity pattern: Calculate on the same loop than statesmap
    StorageIndex const sparcitypattern =
      aVirtualSys.events_.count() * aVirtualSys.states_number_;

    // Map states to its new index
    // TODO: It SHOULD be returned in the future in a pair
    StorageIndex cst = 0;
    for (StorageIndex s : aVirtualSys.virtual_states_) {
        statesmap[s] = cst;
        aVirtualSys.bittriplet_.push_back(BitTriplet(cst, cst, true));
        ++cst;
    }

    // virtual_states_ is not necessary anymore and it can be a large vector
    aVirtualSys.virtual_states_.clear();

    // Remap marked states
    for (StorageIndex s0 : aVirtualSys.sys0_.GetMarkedStates()) {
        for (StorageIndex s1 : aVirtualSys.sys1_.GetMarkedStates()) {
            StorageIndex const key = s1 * aVirtualSys.n_states_sys0_ + s0;
            if (statesmap.contains(key)) {
                aVirtualSys.InsertMarkedState(statesmap[key]);
            }
        }
    }

    // Reserve space for transitions
    aVirtualSys.triplet_.reserve(sparcitypattern);
    aVirtualSys.bittriplet_.reserve(sparcitypattern +
                                    aVirtualSys.states_number_);

    // Calculate transitions
    while (!aVirtualSys.transtriplet_.empty()) {
        auto q_trans = aVirtualSys.transtriplet_.back();
        auto const q = q_trans.first;

        while (!q_trans.second->empty()) {
            auto const qto_e = q_trans.second->back();
            auto const qto = qto_e.first;

            if (statesmap.contains(qto)) {
                auto const event = qto_e.second;
                auto const qto_mapped = statesmap[qto];
                auto const q_mapped = statesmap[q];

                aVirtualSys.triplet_.push_back(Triplet<NEvents>(
                  q_mapped, qto_mapped, EventsSet<NEvents>{ 1ul << event }));
                aVirtualSys.bittriplet_.push_back(
                  BitTriplet(qto_mapped, q_mapped, true));
            }
            q_trans.second->pop_back();
        }
        delete q_trans.second;
        aVirtualSys.transtriplet_.pop_back();
    }

    return;
}

/*! \brief Remove bad states recursively
 * \details Remove a state and all the states that become a bad state when the
 * previous one is removed.
 *
 * @param[out] aVirtualSys reference to the virtual system which will have a
 * state removed
 * @param aC A hash table containing the states currently added to the
 * virtual system.
 * @param aQ The state to be removed.
 * @param aNonContrBit Bit array containing the virtual system non controllable
 * events.
 * @param aRmTable A hash table containing all the removed states so far
 * \return void
 */
template<uint8_t NEvents, typename StorageIndex>
void
RemoveBadStates(SyncSysProxy<NEvents, StorageIndex>& aVirtualSys,
                TransMap<StorageIndex>& aC,
                StorageIndex const& aQ,
                EventsSet<NEvents> const& aNonContrBit,
                StatesTableHost<StorageIndex>& aRmTable)
{
    StatesStack<StorageIndex> f;
    f.push(aQ);
    aRmTable.insert(aQ);

    while (!f.empty()) {
        auto const x = f.top();
        f.pop();

        auto q_events = aVirtualSys.GetInvStateEvents(x);

        q_events &= aNonContrBit;

        cldes::ScalarType event = 0;
        while (q_events.any()) {
            if (q_events.test(0)) {
                StatesArray<StorageIndex> const finv =
                  aVirtualSys.InvTrans(x, event);

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

/*! \brief Computes the monolithic supervisor of a plant and a spec
 *
 * @param aP Plant system const reference
 * @param aE Specs system const reference
 * @param aNonContr Hash table containing all non-controllable events indexes.
 * \return The monolithic supervisorconcrete system
 */
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(DESystemBase<NEvents, StorageIndex> const& aP,
                DESystemBase<NEvents, StorageIndex> const& aE,
                EventsTableHost const& aNonContr)
{
    // Define new systems params: Stage1 is not necessary
    SyncSysProxy<NEvents, StorageIndex> virtualsys{ aP, aE };

    // non_contr in a bitarray structure
    EventsSet<NEvents> non_contr_bit;
    EventsSet<NEvents> p_non_contr_bit;

    // Evaluate which non contr event is in system and convert it to a
    // bitarray
    for (cldes::ScalarType event : aNonContr) {
        if (aP.GetEvents().test(event)) {
            p_non_contr_bit.set(event);
            if (virtualsys.events_.test(event)) {
                non_contr_bit.set(event);
            }
        }
    }

    // Supervisor states
    TransMap<StorageIndex> c;
    StatesTableHost<StorageIndex> rmtable;

    // f is a stack of states accessed in a dfs
    StatesStack<StorageIndex> f;

    // Initialize f and ftable with the initial state
    f.push(virtualsys.init_state_);

    // Allocate inverted graph, since we are search for inverse transitions
    virtualsys.AllocateInvertedGraph();

    while (!f.empty()) {
        auto const q = f.top();
        f.pop();

        if (!rmtable.contains(q) && !c.contains(q)) {
            auto const qx = q % virtualsys.n_states_sys0_;
            auto const q_events = virtualsys.GetStateEvents(q);

            auto const in_ncqx = p_non_contr_bit & aP.GetStateEvents(qx);
            auto const in_ncqx_and_q = in_ncqx & q_events;

            if (in_ncqx_and_q != in_ncqx) {
                // TODO: Fix template implicit instantiation
                RemoveBadStates<NEvents, StorageIndex>(
                  virtualsys, c, q, non_contr_bit, rmtable);
            } else {
                c[q] = new InvArgTrans<StorageIndex>();

                cldes::ScalarType event = 0;
                auto event_it = q_events;
                while (event_it.any()) {
                    if (event_it.test(0)) {
                        auto const fsqe = virtualsys.Trans(q, event);

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

    rmtable.clear();

    virtualsys.ClearInvertedGraph();

    // Swap new system states and sort it
    virtualsys.virtual_states_.reserve(c.size());
    virtualsys.transtriplet_.reserve(c.size());
    for (auto tr : c) {
        virtualsys.virtual_states_.push_back(tr.first);
        virtualsys.transtriplet_.push_back(tr);
    }
    c.clear();

    // Finish synching without removed states
    // SynchronizeStage2(virtualsys);

    // Transform virtual sys in a real sys by forcing conversion
    auto sys = DESystem<NEvents, StorageIndex>(virtualsys);

    // Remove non-accessible and non-coaccessible states
    // TODO: Trim will be performed by GPU
    // sys.Trim();

    // bye
    return sys;
}

/*! \brief Polymorphic shared ptr to base
 */
template<uint8_t NEvents, typename StorageIndex>
using Node = std::shared_ptr<DESystemBase<NEvents, StorageIndex>>;

/*! \brief Vector of polymorphic smart pointers
 */
template<uint8_t NEvents, typename StorageIndex>
using RefNodes = std::vector<Node<NEvents, StorageIndex>>;

/*! \brief Binary tree root and reference to nodes
 */
template<uint8_t NEvents, typename StorageIndex>
using BinExprTree =
  std::pair<Node<NEvents, StorageIndex>, RefNodes<NEvents, StorageIndex>>;

/*! \brief Generate a expression tree of synchronize operations
 * \details Build a binary expression tree of parallel compositions
 * It is not necessary to keep a traditional vector representing
 * a binary tree, since each object keep track of its sons.
 * But we still need to keep track of them, since the SyncSysProxy uses
 * references, not copies. The tree is balanced like the following example:
 * e.g. to a vector of plants of size 5:
 *         P
 *        /\
 *       || p4
 *      / \
 *     /   \
 *    /     \
 *   ||     ||
 *  /  \    / \
 * p0  p1  p2 p3
 *
 * \param[in] aSystems Vector of systems to be mutually synchronized.
 * \return Binary tree pair: root and references
 */
template<uint8_t NEvents, typename StorageIndex>
BinExprTree<NEvents, StorageIndex>
GenBinExprTree(DESVector<NEvents, StorageIndex> const& aSystems)
{
    using SyncSysProxy = SyncSysProxy<NEvents, StorageIndex>;
    using DESystemBase = DESystemBase<NEvents, StorageIndex>;
    using DESystem = DESystem<NEvents, StorageIndex>;

    std::vector<std::shared_ptr<DESystemBase>> sys;
    std::vector<std::shared_ptr<DESystemBase>> nodes_ref;

    for (auto s : aSystems) { // initialize tree
        auto node = std::make_shared<DESystem>(s);
        sys.push_back(node);
        nodes_ref.push_back(node);
    }
    while (sys.size() != 1) {
        auto cp_sys = std::move(sys);
        if (cp_sys.size() % 2 != 0) {
            std::shared_ptr<DESystemBase> node =
              cp_sys.back(); // node is a shared ptr
            cp_sys.pop_back();
            sys.push_back(node);
        }
        size_t processed_nodes = 0ul;
        while (!cp_sys.empty()) { // So it has even number of items
            std::shared_ptr<DESystemBase> lhs = cp_sys.back();
            cp_sys.pop_back();
            std::shared_ptr<DESystemBase> rhs = cp_sys.back();
            cp_sys.pop_back();
            std::shared_ptr<DESystemBase> node =
              std::make_shared<SyncSysProxy>(SyncSysProxy{ *lhs, *rhs });
            nodes_ref.push_back(node);
            sys.push_back(node);
            processed_nodes += 2;
        }
    }
    return std::make_pair(sys[0], nodes_ref);
}

/*! \brief Computes the monolithic supervisor of plants and specs
 * \details Build a binary expression tree of synchronizations and execute
 * the supervisor synthesis with the second level of the tree.
 * The root of the tree is the supervisor.
 * \warning Assumes that vectors has length >= 2
 * \note It balance the trees with a naive algorithm. The most
 * effiecient this structure is, the most scattered the events are
 * on the leafs: Not only the balance matters here.
 *
 * @param aPlants Vector containing plants systems const references
 * @param aSpecs Vector containing specs systems const references
 * @param aNonContr Hash table containing all non-controllable events indexes
 * \return The monolithic supervisor concrete system
 */
template<uint8_t NEvents, typename StorageIndex>
DESystem<NEvents, StorageIndex>
SupervisorSynth(DESVector<NEvents, StorageIndex> const& aPlants,
                DESVector<NEvents, StorageIndex> const& aSpecs,
                EventsTableHost const& aNonContr)
{
    using BinExprTree = BinExprTree<NEvents, StorageIndex>;

    BinExprTree plant = GenBinExprTree(aPlants);
    BinExprTree spec = GenBinExprTree(aSpecs);

    auto const supervisor = SupervisorSynth<NEvents, StorageIndex>(
      *(plant.first), *(spec.first), aNonContr);

    return supervisor;
}

} // namespace op
} // namespace cldes

#endif // DESYSTEM_HPP
