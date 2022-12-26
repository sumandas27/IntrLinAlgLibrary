#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {
    
    auto start = std::chrono::steady_clock::now();

    /* TEST CODE */

    using namespace ila;

    Matrix<4, 4> A
    (
        -3, 2, 0, 0,
        -3, 4, 0, 0,
        0, 0, -5, -4,
        0, 0, -2, 2
    );

    print(generate_eigenvalues(A));

    return 0;
}