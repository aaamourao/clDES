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

 File: cldes/src/des/DESystemCLFwd.hpp
 Description: DESystemCL forward declarations, includes, aliases...
 =========================================================================
*/


#include <CL/cl.hpp>

#define VIENNACL_WITH_EIGEN 1

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/ocl/backend.hpp"

// OclBackend should always be included after viennacl headers
#include "cldes/backend/OclBackend.hpp"

template<typename T>
struct Eigen_vector
{
    typedef typename T::ERROR_NO_EIGEN_TYPE_AVAILABLE error_type;
};
template<>
struct Eigen_vector<float>
{
    typedef Eigen::VectorXf type;
};
template<>
struct Eigen_vector<double>
{
    typedef Eigen::VectorXd type;
};

namespace cldes {

template<uint8_t NEvents, typename StorageIndex>
class DESystem;

namespace op {
namespace backend {
class OclBackend;
}
}
}
