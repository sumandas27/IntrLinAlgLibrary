#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 3> mat1(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0
    );

    Matrix<3, 2> mat2(
        7.0, 10.0,
        8.0, 11.0,
        9.0, 12.0
    );

    Matrix<2, 2> product = mat1 * mat2;

    product.print();

    return 0;
}