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

 File: OclBackend.cpp
 Description: Singleton class implementation: encapsulates the custom
 kernels used by clDES.
 =========================================================================
*/

#include "utils/oclbackend.hpp"

using namespace cldes::utils;

// Initialize static data members
OclBackend* OclBackend::instance_ = nullptr;
viennacl::ocl::program* OclBackend::cldes_program_ = nullptr;
viennacl::ocl::kernel* OclBackend::add_devkernel_ = nullptr;

OclBackend* OclBackend::Instance() {
    if (!instance_) {
        instance_ = new OclBackend();
    }
    return instance_;
}

OclBackend::OclBackend() {
    // All custom clDES OpenCL kernels
    // TODO: put it on other file
    // TODO: This is a dense matrix kernel
    const char* cldes_kernels_ =
        "__kernel void elementwise_add(\n"
        "          __global const float * vec1,\n"
        "          __global const float * vec2, \n"
        "          __global float * result,\n"
        "          unsigned int size) \n"
        "{ \n"
        "  for (unsigned int i = get_global_id(0); i < size; i += "
        "      get_global_size(0))\n"
        "    result[i] = vec1[i] + vec2[i];\n"
        "};\n";

    // Load kernels on device
    cldes_program_ = &viennacl::ocl::current_context().add_program(
        cldes_kernels_, "cldes_kernels_");
}

viennacl::ocl::kernel& OclBackend::AddKernel() {
    if (add_devkernel_) {
        add_devkernel_ = &cldes_program_->get_kernel("elementwise_add");
    }

    return *add_devkernel_;
}

void OclBackend::Enqueue(viennacl::ocl::kernel& aCustomKernel) {
    viennacl::ocl::enqueue(aCustomKernel);
}
