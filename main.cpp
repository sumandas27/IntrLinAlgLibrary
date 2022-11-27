#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    std::array<double, 5> vec1{};
    std::cout << vec1[1];

   /*std::vector<double> matrix
    {
        1, 2, -1, 2, 1, 2,
        -1, -2, 1, 2, 3, 6,
        2, 4, -3, 2, 0, 3,
        -3, -6, 2, 0, 3, 9
    };

    Matrix m = Matrix(4, 6, matrix);
    ERO_row_sum(m, 1, 1, 2);
    ERO_row_sum(m, -2, 1, 3);
    ERO_row_sum(m, 3, 1, 4);
    ERO_row_swap(m, 2, 3);
    ERO_row_sum(m, -1, 2, 4);
    ERO_row_sum(m, -2, 3, 4);
    ERO_scalar_multiplication(m, 1.0/4, 3);
    ERO_row_sum(m, 2, 3, 2);
    ERO_row_sum(m, -2, 3, 1);
    ERO_scalar_multiplication(m, -1, 2);
    ERO_row_sum(m, 1, 2, 1);
    m.print();*/

    return 0;
}