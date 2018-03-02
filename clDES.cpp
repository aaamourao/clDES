// ViennaCL includes
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
// Additional helper functions for this tutorial:
#include "vector-io.hpp"
// Shortcut for writing 'ublas::' instead of 'boost::numeric::ublas::'
using namespace boost::numeric;
int main() {
    typedef float ScalarType;
    std::size_t size = 5;
    ublas::vector<ScalarType> rhs =
        ublas::scalar_vector<ScalarType>(size, ScalarType(size));
    ublas::compressed_matrix<ScalarType> ublas_matrix(size, size);
    ublas_matrix(0, 0) = 2.0f;
    ublas_matrix(0, 1) = -1.0f;
    ublas_matrix(1, 0) = -1.0f;
    ublas_matrix(1, 1) = 2.0f;
    ublas_matrix(1, 2) = -1.0f;
    ublas_matrix(2, 1) = -1.0f;
    ublas_matrix(2, 2) = 2.0f;
    ublas_matrix(2, 3) = -1.0f;
    ublas_matrix(3, 2) = -1.0f;
    ublas_matrix(3, 3) = 2.0f;
    ublas_matrix(3, 4) = -1.0f;
    ublas_matrix(4, 3) = -1.0f;
    ublas_matrix(4, 4) = 2.0f;
    std::cout << "ublas matrix: " << ublas_matrix << std::endl;
    viennacl::vector<ScalarType> vcl_rhs(size);
    viennacl::compressed_matrix<ScalarType> vcl_compressed_matrix(size, size);
    viennacl::copy(rhs, vcl_rhs);
    viennacl::copy(ublas_matrix, vcl_compressed_matrix);
    // just get the data directly from the GPU and print it:
    ublas::compressed_matrix<ScalarType> temp(size, size);
    viennacl::copy(vcl_compressed_matrix, temp);
    std::cout << "ViennaCL: " << temp << std::endl;
    // now modify GPU data directly:
    std::cout << "Modifying vcl_compressed_matrix a bit: " << std::endl;
    vcl_compressed_matrix(0, 0) = 3.0f;
    vcl_compressed_matrix(2, 3) = -3.0f;
    vcl_compressed_matrix(4, 2) = -3.0f; // this is a new nonzero entry
    vcl_compressed_matrix(4, 3) = -3.0f;
    // and print it again:
    viennacl::copy(vcl_compressed_matrix, temp);
    std::cout << "ViennaCL matrix copied to uBLAS matrix: " << temp
              << std::endl;
    std::cout << "ublas: " << ublas::prod(temp, rhs) << std::endl;
    std::cout << "ViennaCL: "
              << viennacl::linalg::prod(vcl_compressed_matrix, vcl_rhs)
              << std::endl;
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}
