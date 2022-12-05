#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    /*Vector<5> v1(1, 2, 3, 4, 5);
    Vector<5> v2(6, 4, 2, 7, 5);
    Vector<5> sum = v1 + v2;

    sum.print();*/

    Vector<4> vec(1, 2, 3, 4);
    Vector<4> zero = zero_vector<4>();
    Vector<4> test = vec + zero;
    test.print();

    return 0;
}