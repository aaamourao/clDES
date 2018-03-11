# clDES

## Introduction

C++11/OpenCL library for **Automata** as a **Discrete-Event Systems** view
and its operations on a **parallel approach** . A **Discrete Controller
Synthesizer** parallel algorithm is also available.

## Available operations

* Accessible Part
* Coaccessible part
* Trim
* Controller Synthesis

## Compiling

It is necessary to make sure that CMake is installed, such as gnu C/C++
developer tools.

### Build project

You gonna need have system support to OpenCL version `>= 1.2`, and Boost
Library version `>= 1.64`. After checking that, compile the project:

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
