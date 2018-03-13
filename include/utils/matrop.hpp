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

 File: matrop.hpp
 Description: OpenCL Sparse Matrix operations not implemented by
 ViennaCL.
 =========================================================================
*/

#ifndef MATRIX_OPERATIONS_HPP
#define MATRIX_OPERATIONS_HPP

#ifndef VIENNACL_WITH_OPENCL
#define VIENNACL_WITH_OPENCL
#endif

#include "constants.hpp"
#include "viennacl/compressed_matrix.hpp"

namespace cldes::utils {

using ScalarType = cldes::ScalarType;

/*! \brief utils::Sum function
 *
 * Sum two sparse matrix located on device memory.
 */
viennacl::compressed_matrix<ScalarType> Sum(
    viennacl::compressed_matrix<ScalarType> &aMatrix1,
    viennacl::compressed_matrix<ScalarType> &aMatrix2);

}  // namespace cldes::utils

#endif  // MATRIX_OPERATIONS_HPP
