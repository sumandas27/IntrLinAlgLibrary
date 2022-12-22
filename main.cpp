#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using namespace ila;

    Matrix<3, 5> m
    (
        -3, 6, -1, 1, -7,
        1, -2, 2, 3, -1,
        2, -4, 5, 8, -4
    );

    std::cout << dim(row(m));

    return 0;
}