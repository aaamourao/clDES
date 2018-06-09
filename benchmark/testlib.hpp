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

 File: testlib.hpp
 Description: Library for common functions used on clDES tests.
 =========================================================================
*/

#include <Eigen/Sparse>
#include <iostream>

template<typename T, typename StringType>
std::string
ReadResult(T const& aOpResult, StringType const aHeader)
{
    std::ostringstream result;

    result << aHeader << ": ";
    for (auto state : aOpResult) {
        result << state << " ";
    }
    result << ">" << std::endl;

    return result.str();
}

using eigen_matrix = Eigen::SparseMatrix<std::bitset<56>, Eigen::RowMajor>;

template<typename StringType>
std::string
ReadResult(eigen_matrix const& aOpResult, StringType const aHeader)
{
    std::ostringstream result;

    result << aHeader << ":" << std::endl;
    for (auto it1 = 0l; it1 != aOpResult.rows(); ++it1) {
        for (auto it2 = 0l; it2 != aOpResult.cols(); ++it2) {
            result << aOpResult.coeff(it1, it2).to_ullong() << " ";
        }
        result << std::endl;
    }
    result << ">" << std::endl;

    return result.str();
}

template<typename StringType>
void
TestResult(StringType const& aResult, StringType const& aExpected)
{
    assert(aResult.find(aExpected) != std::string::npos);
}

template<typename T, typename StringType>
void
ProcessResult(T const& aOpResult,
              StringType const aHeader,
              StringType const aExpected)
{
    auto result = ReadResult(aOpResult, aHeader);
    std::cout << "Result:   " << result;
    std::cout << "Expected: " << aHeader << ": " << aExpected << std::endl;
    TestResult(result, std::string(aExpected));
}

template<typename GraphType, typename StringType>
void
PrintGraph(GraphType const& aGraph, StringType const& aGraphName)
{
    std::cout << aGraphName << std::endl;
    for (auto it1 = 0l; it1 != aGraph.rows(); ++it1) {
        for (auto it2 = 0l; it2 != aGraph.cols(); ++it2) {
            std::cout << aGraph.coeff(it1, it2).to_ullong() << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
