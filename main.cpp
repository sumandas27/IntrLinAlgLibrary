#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    std::array<double, 6> m =
    {
        1.0, 4.0,
        2.0, 5.0,
        3.0, 6.0
    };

    std::array<double, 2> v = { 7.0, 8.0 };

    Matrix<3, 2> mat
    ({
        1.0, 4.0,
        2.0, 5.0,
        3.0, 6.0
    });

    Vector<2> vec({ 7.0, 8.0 });

    Vector<3> product = mat * vec;

    product.print();

    return 0;
}