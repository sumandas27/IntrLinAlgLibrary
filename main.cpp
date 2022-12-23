#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using namespace ila;

    Vector<2> w(4, 4);
    normalize(w);
    print(w);

    return 0;
}