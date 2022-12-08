#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<2, 2> mat(
        3, 2,
        -1, 0
    );

    solve_homogenous_system(5 * mat - 5 * identity_matrix<2>()).print();

    Vector<3> vec(-1, 2, 1);



    return 0;
}