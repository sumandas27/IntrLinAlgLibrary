#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<3, 3> A
    (
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    );

    std::cout << det(A);

    return 0;
}