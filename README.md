# clDES

## Introduction

C++11/OpenCL library for **Automata** as a **Discrete-Event Systems** view
and its operations on a **parallel approach** . A **Discrete Controller
Synthesizer** parallel algorithm is also available.

## Implementation

**clDES** uses a graph based model to represent discrete-event systems,
which are implemented with adjacency sparse matrices. It has a
[ViennaCL](http://viennacl.sourceforge.net/) back-end for multiplying
sparse matrices on GPU and a custom **OpenCL** kernel for the parallel
composition operation.

### Available operations

The following operations are going to be available on **clDES-1.0.0**:

* Accessible part
* Coaccessible part
* Trim
* Automata synchronization: parallel composition (*in progress*)
* Controller Synthesis (*will be implemented soon*)

## Compiling

Follow the instructions bellow to build the project.

### Requirements

It is necessary to make sure that **CMake** is installed, such as **GNU** C/C++
developer tools.

You gonna need have system support to **OpenCL** version `>= 1.2`, and **Boost**
Library version `>= 1.64`. So far, you should install it manually. In a near
future, **CMake** will install it locally if it was not previously installed in
the system.

### Build project

After checking if all requirements are satisfied, compile the project:

```bash
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake ..
$ make
```

### Dev build

If you wanna contribute **clDES**, you may want to compile with debug flags
setting **CMake**'s flag `BUILD_TYPE` to `Debug`:

```bash
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
