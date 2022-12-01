#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    Matrix<3, 4> mat
    ({
        1, 3, 0, 2,
        0, 0, 1, 4,
        1, 3, 1, 6
    });

    Vector<3> vec({1, 6, 5});

    solve(mat, vec).print();
    std::cout << is_consistent(mat, vec);

    return 0;
}