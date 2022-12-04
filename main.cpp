#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<3, 4> mat
    (
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12, 13
    );

    mat.print();

    return 0;
}