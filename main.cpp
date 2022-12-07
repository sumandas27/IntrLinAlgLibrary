#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 2> mat(
        1, 2,
        1, 2
    );

    std::cout << inverse(mat);

    return 0;
}