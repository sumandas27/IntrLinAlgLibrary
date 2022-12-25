#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using namespace ila;

    Matrix<3, 2> m
    (
        1, 2,
        1, 2,
        1, 2
    );

    auto [q, r] = qr_factorization(m);
    ila::print(q);
    std::cout << "\n";
    ila::print(r);

    return 0;
}