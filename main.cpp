#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    /*// test vector
    std::vector<double> testVector1{ 1.0, 2.0, 3.0 , 4.0 };
    std::vector<double> testVector2{ 1.0, 3.0, 10.0, 8.0 };
    Vector v1 = Vector(4, testVector1);
    Vector v2 = Vector(4, testVector2);
    Vector p1 = 2 * v1;
    Vector p2 = v2 * 3;
    Vector sum = p1 + p2;
    sum.print();
    p1.print();
    p2.print();*/

    /*std::vector<double> testMatrix1
    {
        1.0, 2.0, 3.0, 5.0,  20.0,
        4.0, 5.0, 6.0, 10.0, 25.0,
        7.0, 8.0, 9.0, 15.0, 0.0
    };

    std::vector<double> testMatrix2
    {
        1.0, 2.0, 3.0, 5.0,  20.0,
        4.0, 5.0, 6.0, 10.0, 25.0,
        7.0, 8.0, 9.0, 15.0, 0.0
    };

    Matrix m1 = Matrix(3, 5, testMatrix1);
    Matrix m2 = Matrix(3, 5, testMatrix2);
    //m1.set(3, 4, 100.0);
    std::cout << (m1 == m2);*/

    Vector test = zero_vector(5);
    test.set(2, 1.5364);
    test.print();
    zero_matrix(8, 6).print();

    return 0;
}