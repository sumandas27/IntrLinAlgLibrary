#include "IntrLinAlgLibrary.hpp"

int main(int argc, char** argv) {
    
    IntrLinAlgLibrary_init();

    std::cout << EM_scalar_multiplication<5>(3.5, 3);

    return 0;
}