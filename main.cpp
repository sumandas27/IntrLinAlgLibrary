#include "IntrLinAlgLibrary.hpp"

//TODO: Test matrix == operator overload

int main(int argc, char** argv) {
    // test vector
    std::vector<double> testVector1{ 1.0, 2.0, 3.0, 4.0, 7.0 };
    std::vector<double> testVector2{ 1.0, 2.0, 3.0, 4.0, 7.0 };
    Vector v1 = Vector(5, testVector1);
    Vector v2 = Vector(5, testVector2);
    //v1.set(1, 2.0);
    std::cout << (v1 == v2) << "\n";

    // test matrix
    std::vector< std::vector<double> > testMatrix1
    {
        { 1.0, 2.0, 3.0, 5.0,  20.0 },
        { 4.0, 5.0, 6.0, 10.0, 25.0 },
        { 7.0, 8.0, 9.0, 15.0, 0.0  }
    };

    std::vector< std::vector<double> > testMatrix2
    {
        { 1.0, 2.0, 3.0, 5.0,  20.0 },
        { 4.0, 5.0, 6.0, 10.0, 25.0 },
        { 7.0, 8.0, 9.0, 15.0, 0.0  }
    };
    
    Matrix m1 = Matrix(3, 5, testMatrix1);
    Matrix m2 = Matrix(3, 5, testMatrix2);
    //m2.set(2, 2, 3.0);
    std::cout << (m1 == m2) << "\n";
    return 0;
}