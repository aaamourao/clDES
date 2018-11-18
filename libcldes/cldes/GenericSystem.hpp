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

#include <typeinfo>

template<uint8_t NEvents, typename StorageIndex, class RealDESystem>
class DESystemBase;

namespace cldes {
struct GenericSystem
{
    struct InnerSystemBase
    {
        using ptr = std::unique_ptr<InnerSystemBase>;

        virtual ~InnerSystemBase() = default;

        virtual std::type_info const& type() const noexcept = 0;

        virtual InnerSystemBase* clone() const noexcept = 0;

        virtual long unsigned size() const noexcept = 0;

        virtual bool isVirtual() const noexcept = 0;
    };

    template<typename SysT_>
    struct InnerSystem : InnerSystemBase
    {
        /* User should define:
         * EventsSet_t
         * StatesSet_t
         * Event_t
         * State_t
         * Size_t
         */

        InnerSystem(SysT_ aSys)
          : innersys_{ std::move(aSys) }
        {}

        virtual std::type_info const& type() const noexcept override
        {
            return typeid(SysT_);
        }

        virtual InnerSystemBase* clone() const noexcept override
        {
            return new InnerSystem{ SysT_{ innersys_ } };
        }

        SysT_& operator*() { return innersys_; }

        SysT_ const& operator*() const { return innersys_; }

        virtual long unsigned size() const noexcept override;

        virtual bool isVirtual() const noexcept override;

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

    long unsigned size() const noexcept { return inner_->size(); }

    bool isVirtual() const noexcept { return inner_->isVirtual(); }

    typename InnerSystemBase::ptr inner_;
};

template<class SysT_l, class SysT_r>
bool constexpr
operator==(GenericSystem const& aLhs, GenericSystem const& aRhs)
{
    return aLhs.template cast<SysT_l>() == aRhs.template cast<SysT_r>();
}

template<>
struct SysTraits<GenericSystem>
{
    uint8_t static constexpr Ne_ = 64;
    using Si_ = unsigned long;
};
}
#endif // GENERIC_SYSTEM_HPP
