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

#include "backend/matrop.hpp"
#include "backend/oclbackend.hpp"

using namespace cldes::backend;

viennacl::compressed_matrix<ScalarType> Sum(
    viennacl::compressed_matrix<ScalarType> &aMatrix1,
    viennacl::compressed_matrix<ScalarType> &aMatrix2) {
    auto oclbackend = OclBackend::Instance();
    auto &my_kernel_mul = oclbackend->AddKernel();

    viennacl::compressed_matrix<ScalarType> result_sum(aMatrix1.size1(),
                                                       aMatrix1.size2());
    oclbackend->Enqueue(
        my_kernel_mul(aMatrix1, aMatrix2, result_sum,
                      static_cast<cl_uint>(aMatrix1.size1())));

    return result_sum;
}
