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

#include "backend/oclbackend.hpp"
#include <CL/cl.hpp>
#include <fstream>
#include <string>

using namespace cldes::backend;

// Initialize static data members
OclBackend* OclBackend::instance_ = nullptr;

OclBackend* OclBackend::Instance() {
    if (!instance_) {
        instance_ = new OclBackend();
    }
    return instance_;
}

OclBackend::OclBackend() {
    context_ = viennacl::ocl::current_context();

    // Read kernels from file
    std::ifstream source_file_name("kernels.cl");

    std::string source_file(std::istreambuf_iterator<char>(source_file_name),
                            (std::istreambuf_iterator<char>()));

    // Load all custom cldes kernels on device
    cldes_program_ = context_.add_program(source_file, "cldes_kernels");
}

OclBackend::ViennaCLContext OclBackend::GetContext() {
    return context_;
}

OclBackend::ViennaCLKernel& OclBackend::GetKernel(std::string aKernelName) {
    return cldes_program_.get_kernel(aKernelName);
}

void OclBackend::Enqueue(viennacl::ocl::kernel const& k) {
    viennacl::ocl::enqueue(k);
    queue_ = viennacl::ocl::get_queue().handle().get();
}

OclBackend::CLQueue OclBackend::CommandQueue() { return queue_; }
