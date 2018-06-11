#clDES

##Introduction

  C++ 11/*OpenCL 1.2* library for **Discrete-Event Systems** modeled as **Automata**
and their operations on a **parallel programming approach** . A
**Supervisor Synthesizer** parallel algorithm is also available.

## Implementation

**clDES** uses a graph based model to represent discrete-event systems,
which are implemented with adjacency sparse matrices. It has a
[ViennaCL](http://viennacl.sourceforge.net/) back-end for multiplying
sparse matrices on GPU and customs **OpenCL** kernels for operations
such as the parallel composition.

### Available operations

> Only CPU operations are available now. The *OpenCL* code is being refactored.

The following operations are going to be available on **clDES-1.0.0**:

Automata operation | Implementation | Status
-------------------|----------------|--------
Accessible part | `cldes::DESystem<NEvents, StorageIndex>::AccessiblePart()` | **Implemented**
Coaccessible part | `cldes::DESystem<NEvents, StorageIndex>::CoaccessiblePart()` | **Implemented**
Trim | `cldes::DESystem<NEvents, StorageIndex>::Trim()` | **Implemented**
Synchronization: parallel composition | `cldes::op::Synchronize<NEvents, StorageIndex>()` | **Implemented**
Virtual parallel composition | `cldes::op::SynchronizeStage1<NEvents, StorageIndex>()` | **Implemented**
Real parallel composition | `cldes::op::SynchronizeStage2<NEvents, StorageIndex>()` | **Implemented**
Controller Synthesis | `cldes::op::SupervisorSynth<NEvents, StorageIndex>()` | **Implemented**

The following operations are going to be available on **clDES-1.1.0**:

Implementation of graphs representing **DES** as **tensors** instead of
matrices: an even more efficiente parallel algorithm for
`cldes::op::Synchronize` will be possible.

## Compiling

Follow the instructions bellow to build the project.

### Requirements

It is necessary to make sure that **CMake 3.11** is installed, such as
**LLVM Clang 5.0.2** or **GCC/G++ 7.3.1** developer tools.

and **Eigen3** Library
are required. So far, you should install it manually. In a near
future, **CMake** will install it locally if it was not previously installed in
the system.

> It may compile and run smoothly on older versions. These tools mentioned above
> compose the tool chain used to develop and to test so far (**Fedora 27**).

System support to **OpenCL** version `>= 1.2`, **Boost** Library version `>= 1.58`
**clDES** embedded libraries so far:

* [ViennaCL 1.7.1](Efficient Linear Algebra Algorithms)
* [Sparsepp](https://github.com/greg7mdp/sparsepp): The faster sparse hash
  sets/maps **C++11** implementation I know.

> It is always necessary to set `CLDES_INCLUDE_PATH` to run **clDES** applications
> which execute operations that use custom *OpenCL* kernels, such as
> `cldes::op::Syncronization()`.

### Build project

After checking if all requirements are satisfied, compile the project:

```bash
$ export CC=/usr/bin/clang
$ export CXX=/usr/bin/clang++
$ export CLDES_INCLUDE_PATH=<clDES_root>/include/
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake ..
$ make
```

### Dev build

If you wanna contribute to **clDES**, you may want a debug build by setting
**CMake**'s flag `BUILD_TYPE` to `Debug`:

```bash
$ export CC=/usr/bin/clang
$ export CXX=/usr/bin/clang++
$ export CLDES_INCLUDE_PATH=<clDES_root>/include/
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make
```

## Running tests

After inserting code, or if you just want to test the library, execute the tests
using **CTest**:

```
$ cd <clDES_root>/build
$ make test
```
