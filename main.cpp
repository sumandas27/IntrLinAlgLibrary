#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using namespace ila;

    std::array<Vector<4>, 3> basis =
    {
        Vector<4>(1, -1, 0, 2),
        Vector<4>(1, 1, 1, 3),
        Vector<4>(3, 1, 1, 5)
    };
    
    ila::print(basis);
    std::cout << "\n";
    ila::print(orthonormal_basis(basis));

    return 0;
}