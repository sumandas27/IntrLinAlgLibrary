#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    /* TEST CODE */

    using ila::Vector, ila::Matrix;

    std::vector<Vector<3>> mySet =
    {
        Vector<3>(1000, 2, 3),
        Vector<3>(4, 5, 6)
    };

    ila::print(mySet);

    return 0;
}