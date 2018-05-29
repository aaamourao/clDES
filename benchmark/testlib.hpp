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

 Copyright (c) 20ul18 - Adriano Mourao <adrianomourao@protonmail.com>
                      madc0ulww @ [https://github.com/madc0ulww]

 LacSED - Laboratório de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: testlib.hpp
 Description: Library for common functions used on clDES tests.
 =========================================================================
*/

#include <boost/numeric/ublas/matrix_sparse.hpp>

template <typename T, typename StringType>
std::string ReadResult(T const &aOpResult, StringType const aHeader) {
    std::ostringstream result;

    result << aHeader << ": ";
    for (auto state : aOpResult) {
        result << state << " ";
    }
    result << ">" << std::endl;

    return result.str();
}

using ublas_matrix = boost::numeric::ublas::compressed_matrix<std::bitset<100>>;

template <typename StringType>
std::string ReadResult(ublas_matrix const &aOpResult,
                       StringType const aHeader) {
    std::ostringstream result;

    result << aHeader << ":" << std::endl;
    for (auto it1 = 0ul; it1 != aOpResult.size1(); ++it1) {
        for (auto it2 = 0ul; it2 != aOpResult.size2(); ++it2) {
            auto events = aOpResult(it1, it2);
            if (events != 0ul) {
                result << events.to_ulong() << " ";
            } else {
                result << "0 ";
            }
        }
        result << std::endl;
    }
    result << ">" << std::endl;

    return result.str();
}

template <typename StringType>
void TestResult(StringType const &aResult, StringType const &aExpected) {
    assert(aResult.find(aExpected) != std::string::npos);
}

template <typename T, typename StringType>
void ProcessResult(T const &aOpResult, StringType const aHeader,
                   StringType const aExpected) {
    auto result = ReadResult(aOpResult, aHeader);
    std::cout << "Result:   " << result;
    std::cout << "Expected: " << aHeader << ": " << aExpected << std::endl;
    TestResult(result, std::string(aExpected));
}

template <typename GraphType, typename StringType>
void PrintGraph(GraphType const &aGraph, StringType const &aGraphName) {
    std::cout << aGraphName << std::endl;
    for (auto it1 = 0ul; it1 != aGraph.size1(); ++it1) {
        for (auto it2 = 0ul; it2 != aGraph.size2(); ++it2) {
            auto events = aGraph(it1, it2);
            if (events != 0ul) {
                std::cout << events.to_ulong() << " ";
            } else {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
