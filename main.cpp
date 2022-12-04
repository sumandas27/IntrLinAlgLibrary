#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    /*Vector<5> v1(1, 2, 3, 4, 5);
    Vector<5> v2(6, 4, 2, 7, 5);
    Vector<5> sum = v1 + v2;

    sum.print();*/

    Matrix<10, 9> vec;
    vec.print();

    return 0;
}