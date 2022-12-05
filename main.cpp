#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 2> mat(
        0.0, 1.0,
        1.0, 1.0
    );

    std::cout << is_symmetric(mat);

    return 0;
}