#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    VectorSet<5, 4> set
    (
        Vector<5>(1, 2, 3, 4, 5),
        Vector<5>(6, 7, 8, 9, 10),
        Vector<5>(11, 12, 13, 14, 15),
        Vector<5>(16, 17, 18, 19, 20)
    );

    augment_vector_set(set).print();

    return 0;
}