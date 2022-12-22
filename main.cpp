#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using namespace ila;

    Vector<3> v(1, 2, 3);
    v[1] = 4;
    ila::print(v);

    return 0;
}