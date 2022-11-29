#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<4, 6> mat
    ({
        1, 2, -1, 2, 1, 2,
        -1, -2, 1, 2, 3, 6,
        2, 4, -3, 2, 0, 3,
        -3, -6, 2, 0, 3, 9
    });

    ref(mat).print();

    return 0;
}