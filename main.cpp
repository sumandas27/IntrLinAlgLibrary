#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    // test vector
    std::vector<double> testVector{ 1.0, 2.0, 3.0, 4.0, 7.0 };
    Vector v = Vector(5, testVector);
    std::cout << v.get(4) << "\n";

    // test matrix
    std::vector< std::vector<double> > testMatrix
    {
        { 1.0, 2.0, 3.0, 5.0,  20.0 },
        { 4.0, 5.0, 6.0, 10.0, 25.0 },
        { 7.0, 8.0, 9.0, 15.0, 30.0 }
    };
    
    Matrix m = Matrix(3, 5, testMatrix);
    m.set(3, 4, 2.0);
    std::cout << m.get(3, 4) << "\n";
    return 0;
}