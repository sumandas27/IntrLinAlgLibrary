#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 5> a
    (
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10
    );

    a.column_vector(2).print();

    return 0;
}