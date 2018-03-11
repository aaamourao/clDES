# clDES

## Introduction

C++ Implementation of **Automata** and **Discrete-Event Systems**, their operations and **Discrete
Controller Synthesis** on a **parallel approach**.

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

### Development Build

If you wanna contribute **clDES**, you may want to compile with debug flags
setting cmake's flag `BUILD_TYPE` to `Debug`:

```bash
$ cd  <clDES_root>
$ mkdir build; cd build
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make
```
