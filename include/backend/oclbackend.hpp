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

 File: OclBackend.hpp
 Description: Singleton class which encapsulates the custom kernels
 used by clDES.
 =========================================================================
*/

#include "viennacl/ocl/backend.hpp"

namespace cldes::backend {

class OclBackend {
public:
    using ViennaCLKernel = viennacl::ocl::kernel;
    using ViennaCLProgram = viennacl::ocl::program;

    /*! \brief OclBackend::Instance() method
     *
     * Public method that guarantees that only one object of OclBackend is
     * created. Returns a raw pointer to a unique OclBackend.
     */
    static OclBackend* Instance();

    /*! \brief OclBackend::AddKernel() method
     *
     * Public method that returns the kernel for ADD operation between two
     * sparse matrices.
     */
    static viennacl::ocl::kernel& AddKernel();

    /*! \brief OclBackend::Enqueuel() method
     *
     * Public method for executing kernels on the device. It is basically
     * viennacl opencl backend encapsulated.
     */
    static void Enqueue(ViennaCLKernel& aCustomKernel);

protected:
    /*! \brief OclBackend constructor
     *
     * Loads kernels and gets an OpenCL program.
     */
    OclBackend();

private:
    /*! \brief OclBackend::instance_ raw pointer
     *
     * Keeps track of the unique OclBackend object.
     */
    static OclBackend* instance_;

    /*! \brief OclBackend::cldes_program_ raw pointer
     *
     * Kernels loaded on the device.
     */
    static ViennaCLProgram* cldes_program_;

    /*! \brief OclBackend::add_devkernel_ raw pointer
     *
     * Kernel on device memory for add operation.
     */
    static ViennaCLKernel* add_devkernel_;
};

}  // namespace cldes::backend
