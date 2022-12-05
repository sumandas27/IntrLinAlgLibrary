#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<3, 3> mat
    (
        2, 1, 5,
        -1, 0, -1,
        1, 2, 7
    );

    std::cout << solve_homogenous_system(mat);
    return 0;
}