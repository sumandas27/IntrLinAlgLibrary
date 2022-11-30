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

    Vector<5> vec({1.0, 2.0, 3.0, 4.0, -50.0});
    std::cout << vec;

    rref(mat).print();
    std::cout << rref(mat);
    return 0;
}