#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {

    /* TEST CODE */

    using namespace ila;

    ila::Matrix<3, 3> mySquareMat
    (
        5, -10, -5,
        2, 14, 2,
        -4, -8, 6
    );
    
    ila::print(ila::generate_eigenvalues(mySquareMat));

    return 0;
}