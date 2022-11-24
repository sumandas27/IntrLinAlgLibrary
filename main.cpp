#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /*std::vector<double> matrix
    {
        1.0, 4.0,
        2.0, 5.0,
        3.0, 6.0
    };

    std::vector<double> vector{ 1, 0 };

    Matrix a = rotation_matrix(450);
    Vector v = Vector(2, vector);
    Vector product = a * v;*/
    
    Matrix a = rotation_matrix(46);
    Matrix b = rotation_matrix(46);
    std::cout << (a == b);

    return 0;
}