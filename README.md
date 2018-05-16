# clDES

## Introduction

C++11/*OpenCL 1.2* library for **Discrete-Event Systems** modeled as **Automata**
and their operations on a **parallel programming approach** . A
**Discrete Controller Synthesizer** parallel algorithm is also available.

## Implementation

**clDES** uses a graph based model to represent discrete-event systems,
which are implemented with adjacency sparse matrices. It has a
[ViennaCL](http://viennacl.sourceforge.net/) back-end for multiplying
sparse matrices on GPU and customs **OpenCL** kernels for operations
such as the parallel composition.

### Available operations

The following operations are going to be available on **clDES-1.0.0**:

Automata operation | Implementation | Status
-------------------|----------------|--------
Accessible part | `cldes::DESystem::AccessiblePart()` | **Implemented**
Coaccessible part | `cldes::DESystem::CoaccessiblePart()` | **Implemented**
Trim | `cldes::DESystem::Trim()` | **Implemented**
Synchronization: parallel composition | `cldes::op::Synchronize()` | **Implemented**
Virtual parallel composition | `cldes::op::SynchronizeStage1()` | **Implemented**
Real parallel composition | `cldes::op::SynchronizeStage2()` | **Implemented**
Controller Synthesis | `cldes::op::ControllerSynth()` | **in progress**

The following operations are going to be available on **clDES-1.1.0**:

Implementation of graphs representing **DES** as **tensors** instead of
matrices: an even more efficiente parallel algorithm for
`cldes::op::Synchronize` will be possible.

## Compiling

Follow the instructions bellow to build the project.

### Requirements

It is necessary to make sure that **CMake** is installed, such as **LLVM Clang** C/C++
developer tools.

System support to **OpenCL** version `>= 1.2` and **Boost** Library version `>= 1.64`
are required. So far, you should install it manually. In a near
future, **CMake** will install it locally if it was not previously installed in
the system.

### Build project

After checking if all requirements are satisfied, compile the project:

```bash
$ export CC=/usr/bin/clang
$ export CXX=/usr/bin/clang++
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
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make
```

## Running tests

You will need to copy the `*.cl` files from the test directory into the directory where you will run the tests from:

```bash
$ cp <clDES_root>/include/backend/kernels.cl <clDES_root>/build/
$ cp <clDES_root>/include/backend/kernels.cl <clDES_root>/build/bin/tests/
```

On the **clDES** releases, the kernels are going to be loaded from strings
defined on the code. However, on developer versions, such as the version
available now, this step will always be necessary.

After inserting code, or if you just want to test the library, execute the tests
using **CTest**:

```
$ cd <clDES_root>/build
$ make test
```
