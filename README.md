# clDES

## Introduction

C++11/OpenCL library for **Automata** as a **Discrete-Event Systems** view
and its operations on a **parallel approach** . A **Discrete Controller
Synthesizer** parallel algorithm is also available.

## Available operations

* Accessible part
* Coaccessible part
* Trim
* Automata synchronization: parallel composition
* Controller Synthesis

## Compiling

It is necessary to make sure that **CMake** is installed, such as **GNU** C/C++
developer tools.

### Build project

You gonna need have system support to OpenCL version `>= 1.2`, and Boost
Library version `>= 1.64`. So far, you should install it manually. In a near
future, **CMake** will install it locally if it is not previously installed in
the system.

After checking that, compile the project:

```bash
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake ..
$ make
```

### Dev build

If you wanna contribute **clDES**, you may want to compile with debug flags
setting cmake's flag `BUILD_TYPE` to `Debug`:

```bash
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make
```

## Run tests

After inserting code, or if you just want to test the library, run the tests
using **CTest**:

```
$ cd <clDES_root>/build
$ make test
```
