#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 2> A
    (
        1, 2,
        2, 3
    );

    Matrix<2, 3> B
    (
        1, -1, 2,
        1, 0, 1
    );

    Matrix<2, 3> product = inverse(A) * B;
    product.print();

    return 0;
}