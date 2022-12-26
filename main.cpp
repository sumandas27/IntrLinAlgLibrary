#include "IntrLinAlgLibrary.hpp"
#include <chrono>

int main(int argc, char** argv) {

    /* TEST CODE */

    using namespace ila;

    std::vector<int> yo{};
    yo.reserve(4);
    yo.emplace_back(3);
    std::cout << yo.size();

    return 0;
}