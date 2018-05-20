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

#include "constants.hpp"
#include "viennacl/ocl/backend.hpp"

namespace cl {
class CommandQueue;
}

namespace cldes {
namespace backend {

class OclBackend {
public:
    using ViennaCLKernel = viennacl::ocl::kernel;
    using ViennaCLProgram = viennacl::ocl::program;
    using CLQueue = cl_command_queue;
    using ViennaCLContext = viennacl::ocl::context;

    /*! \brief Instantiate clDES OpenCL backend and returns unique object
     *
     * Public method that guarantees that only one object of OclBackend is
     * created. Returns a raw pointer to a unique OclBackend.
     */
    static OclBackend* Instance();

    /*! \brief Return opencl kernel handle
     *
     * Public method that returns the kernel for ADD operation between two
     * sparse matrices.
     *
     * @param aKernelName String containing kernel name
     */
    ViennaCLKernel& GetKernel(std::string aKernelName);

    /*! \brief Execute kernel on device
     *
     * Public method for executing kernels on the device. It is basically
     * viennacl opencl backend encapsulated.
     *
     * @param k OpenCL Kernel
     */
    void Enqueue(viennacl::ocl::kernel const& k);

    /*! \brief Return device OpenCL context
     *
     * Returns the current context created by the ViennaCL backend
     */
    ViennaCLContext GetContext();

    /*! \brief Read buffer and return success code
     *
     *
     * @param k OpenCL Kernel
     */
    CLQueue CommandQueue();

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
    ViennaCLProgram cldes_program_;

    /*!
     *
     */
    CLQueue queue_;

    /*!
     *
     */
    ViennaCLContext context_;
};

}  // namespace backend
}  // namespace cldes
