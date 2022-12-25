#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {
    
    auto start = std::chrono::steady_clock::now();

    /* TEST CODE */

    using namespace ila;

    Matrix<3, 3> m
    (
        5, -10, -5,
        2, 14, 2,
        -4, -8, 6
    );

    std::vector<Eigenvalue> eigenvalues = generate_eigenvalues(m);
    for (const auto& [eigenvalue, multiplicity] : eigenvalues)
        std::cout << eigenvalue << "\tMultiplicity: " << multiplicity << "\n";

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << std::chrono::duration<double, std::milli> (diff).count() << " ms\n";

    return 0;
}