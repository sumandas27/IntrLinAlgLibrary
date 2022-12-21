#include "IntrLinAlgLibrary.hpp"
#include <memory>

void divide() {

    std::unique_ptr<int> yo1 = std::make_unique<int>(10);
    std::unique_ptr<int> yo2 = std::make_unique<int>(5);

    if (*yo1 > 100 || *yo2 > 100)
        return;

    if (*yo1 > 100 || *yo2 > 100)
        return;

    if (*yo1 > 100 || *yo2 > 100)
        return;

    std::cout << *yo1 / *yo2;
}

int main(int argc, char** argv) {
    divide();
    return 0;
}