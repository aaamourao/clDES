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

#include <boost/numeric/ublas/matrix_sparse.hpp>

template <typename T, typename StringType>
std::ostringstream ReadResult(T const &aOpResult, StringType const aHeader) {
    std::ostringstream result;

    result << aHeader << ": ";
    for (auto state : aOpResult) {
        result << state << " ";
    }
    result << std::endl;

    return result;
}

using ublas_matrix = boost::numeric::ublas::compressed_matrix<float>;

template <typename StringType>
std::ostringstream ReadResult(ublas_matrix const &aOpResult,
                              StringType const aHeader) {
    std::ostringstream result;

    auto last_row = 0;
    for (auto it1 = aOpResult.begin1(); it1 != aOpResult.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            if (last_row != it1.index1()) {
                ++last_row;
                result << std::endl;
            }
            result << *it2 << " ";
        }
    }
    result << std::endl;

    return result;
}

template <typename StringType>
void TestResult(StringType const &aResult, StringType const &aExpected) {
    assert(aResult.find(aExpected) != std::string::npos);
}

template <typename T, typename StringType>
void ProcessResult(T const &aOpResult, StringType const aHeader,
                   StringType const aExpected) {
    auto result = ReadResult(aOpResult, aHeader);
    std::cout << result.str();
    std::ostringstream expected;
    TestResult(result.str(), expected.str());
}

template <typename GraphType, typename StringType>
void PrintGraph(GraphType const &aGraph, StringType const &aGraphName) {
    std::cout << aGraphName << std::endl;
    auto last_row = 0;
    for (auto it1 = aGraph.begin1(); it1 != aGraph.end1(); ++it1) {
        for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
            if (last_row != it1.index1()) {
                ++last_row;
                std::cout << std::endl;
            }
            std::cout << *it2 << " ";
        }
    }
    std::cout << std::endl;
}
