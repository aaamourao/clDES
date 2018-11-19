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

 File: cldes/GenericSystem.hpp
 description: Generic wrapper dispatcher to implement a type erasure
 approach, since DESystemBase implement a static polymorphism approach.
 =========================================================================
*/
/*!
 * \file cldes/GenericSystem.hpp
 *
 * \author Adriano Mourao \@madc0ww
 * \date 2018-11-16
 *
 * Generic wrapper dispatcher to implement a type erasure
 * approach, since DESystemBase implement a static polymorphism approach.
 */

#ifndef GENERIC_SYSTEM_HPP
#define GENERIC_SYSTEM_HPP

#include <set>
#include <typeinfo>

namespace cldes {
template<uint8_t NEvents, typename StorageIndex>
struct GenericSystem
{
    using StatesSet = std::set<StorageIndex>;
    using StorageIndexSigned = typename std::make_signed<StorageIndex>::type;
    using StatesArray = std::vector<StorageIndex>;

    struct InnerSystemBase
    {
        using ptr = std::unique_ptr<InnerSystemBase>;

        virtual ~InnerSystemBase() = default;

        virtual std::type_info const& type() const noexcept = 0;

        virtual InnerSystemBase* clone() const noexcept = 0;

        virtual StorageIndex size() const noexcept = 0;

        virtual bool isVirtual() const noexcept = 0;

        virtual StorageIndex getStatesNumber() const noexcept = 0;

        virtual EventsSet<NEvents> getEvents() const noexcept = 0;

        virtual StorageIndex getInitialState() const noexcept = 0;

        virtual StatesSet getMarkedStates() const noexcept = 0;

        virtual EventsSet<NEvents> getStateEvents(StorageIndex const& aQ) const
          noexcept = 0;

        virtual EventsSet<NEvents> getInvStateEvents(
          StorageIndex const& aQ) const noexcept = 0;

        virtual StorageIndexSigned trans(StorageIndex const& aQ,
                                         ScalarType e) const noexcept = 0;

        virtual StatesArray invtrans(StorageIndex const& aQ, ScalarType e) const
          noexcept = 0;

        virtual bool containstrans(StorageIndex const& aQ, ScalarType e) const
          noexcept = 0;

        virtual bool containsinvtrans(StorageIndex const& aQ,
                                      ScalarType e) const noexcept = 0;

        virtual void allocateInvertedGraph() const noexcept = 0;

        virtual void clearInvertedGraph() const noexcept = 0;
    };

    template<typename SysT_>
    struct InnerSystem : InnerSystemBase
    {
        InnerSystem(SysT_ aSys)
          : innersys_{ std::move(aSys) }
        {}

        virtual std::type_info const& type() const noexcept override
        {
            return typeid(SysT_);
        }

        virtual InnerSystemBase* clone() const noexcept override
        {
            return new InnerSystem<SysT_>{ SysT_{ innersys_ } };
        }

        SysT_& operator*() { return innersys_; }

        SysT_ const& operator*() const { return innersys_; }

        StorageIndex size() const noexcept override { return innersys_.size(); }

        bool isVirtual() const noexcept override
        {
            return innersys_.isVirtual();
        }

        StorageIndex getStatesNumber() const noexcept override
        {
            return innersys_.getStatesNumber();
        }

        EventsSet<NEvents> getEvents() const noexcept override
        {
            return innersys_.getEvents();
        }

        StorageIndex getInitialState() const noexcept override
        {
            return innersys_.getInitialState();
        }

        StatesSet getMarkedStates() const noexcept override
        {
            return innersys_.getMarkedStates();
        }

        EventsSet<NEvents> getStateEvents(StorageIndex const& aQ) const
          noexcept override
        {
            return innersys_.getStateEvents(aQ);
        }

        EventsSet<NEvents> getInvStateEvents(StorageIndex const& aQ) const
          noexcept override
        {
            return innersys_.getInvStateEvents(aQ);
        }

        StorageIndexSigned trans(StorageIndex const& aQ, ScalarType e) const
          noexcept override
        {
            return innersys_.trans(aQ, e);
        }

        StatesArray invtrans(StorageIndex const& aQ, ScalarType e) const
          noexcept override
        {
            return innersys_.invtrans(aQ, e);
        }

        bool containstrans(StorageIndex const& aQ, ScalarType e) const
          noexcept override
        {
            return innersys_.containstrans(aQ, e);
        }

        bool containsinvtrans(StorageIndex const& aQ, ScalarType e) const
          noexcept override
        {
            return innersys_.containsinvtrans(aQ, e);
        }

        void allocateInvertedGraph() const noexcept override
        {
            return innersys_.allocateInvertedGraph();
        }

        void clearInvertedGraph() const noexcept override
        {
            return innersys_.clearInvertedGraph();
        }

    private:
        SysT_ innersys_;
    };

    template<class SysT_>
    GenericSystem(SysT_ aSys)
      : inner_{ new InnerSystem<SysT_>{ std::forward<SysT_>(aSys) } }
    {}

    GenericSystem(GenericSystem const& aSys)
      : inner_{ aSys.inner_->clone() }
    {}

    template<typename SysT_>
    GenericSystem& operator=(SysT_ aSys)
    {
        inner_ =
          std::make_unique<InnerSystem<SysT_>>(std::forward<SysT_>(aSys));
        return *this;
    }

    GenericSystem& operator=(GenericSystem const& aSys)
    {
        GenericSystem tmp{ aSys };
        inner_ = std::move(tmp.inner_);
        return *this;
    }

    template<typename SysT_>
    SysT_& cast()
    {
        return *dynamic_cast<InnerSystem<SysT_>&>(*inner_);
    }

    template<typename SysT_>
    SysT_ const& cast() const
    {
        return *dynamic_cast<InnerSystem<SysT_>&>(*inner_);
    }

    std::type_info const& type() const noexcept { return inner_->type(); }

    StorageIndex size() const noexcept { return inner_->size(); }

    bool isVirtual() const noexcept { return inner_->isVirtual(); }

    StorageIndex getStatesNumber() const noexcept
    {
        return inner_->getStatesNumber();
    }

    EventsSet<NEvents> getEvents() const noexcept
    {
        return inner_->getEvents();
    }

    StorageIndex getInitialState() const noexcept
    {
        return inner_->getInitialState();
    }

    StatesSet getMarkedStates() const noexcept
    {
        return inner_->getMarkedStates();
    }

    EventsSet<NEvents> getStateEvents(StorageIndex const& aQ) const noexcept
    {
        return inner_->getStateEvents(aQ);
    }

    EventsSet<NEvents> getInvStateEvents(StorageIndex const& aQ) const noexcept
    {
        return inner_->getInvStateEvents(aQ);
    }
    void allocateInvertedGraph() const noexcept
    {
        return inner_->allocateInvertedGraph();
    }

    StorageIndexSigned trans(StorageIndex const& aQ, ScalarType e) const
      noexcept
    {
        return inner_->trans(aQ, e);
    }

    StatesArray invtrans(StorageIndex const& aQ, ScalarType e) const noexcept
    {
        return inner_->invtrans(aQ, e);
    }

    bool containstrans(StorageIndex const& aQ, ScalarType e) const noexcept
    {
        return inner_->containstrans(aQ, e);
    }

    bool containsinvtrans(StorageIndex const& aQ, ScalarType e) const noexcept
    {
        return inner_->containsinvtrans(aQ, e);
    }

    void clearInvertedGraph() const noexcept
    {
        return inner_->clearInvertedGraph();
    }

    typename InnerSystemBase::ptr inner_;
};

template<class SysT_l,
         class SysT_r,
         uint8_t NEvents = SysTraits<SysT_l>::Ne_,
         typename StorageIndex = typename SysTraits<SysT_l>::Si_>
bool constexpr
operator==(GenericSystem<NEvents, StorageIndex> const& aLhs,
           GenericSystem<NEvents, StorageIndex> const& aRhs)
{
    return aLhs.template cast<SysT_l>() == aRhs.template cast<SysT_r>();
}
}
#endif // GENERIC_SYSTEM_HPP
