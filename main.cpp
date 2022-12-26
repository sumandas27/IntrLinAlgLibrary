#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {

    /* TEST CODE */

    using namespace ila;

    Matrix<3, 3> A
    (
        1, 0, 2,
        0, 3, 0,
        2, 0, 1
    );

    print(get_eigenspace_basis(A, 3));

    return 0;
}