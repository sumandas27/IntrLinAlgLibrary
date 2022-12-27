#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {

    /* TEST CODE */

    using namespace ila;

    Matrix<3, 3> a
    (
        2, 0, 0,
        1, 2, 1,
        -1, 0, 1
    );

    auto [p, d] = diagonalize(a);
    print(p);
    std::cout << "\n";
    print(d);
    std::cout << "\n";
    
    Matrix test = p * d * inverse(p);
    print(test);

    return 0;
}