#include <iostream>
#include "ProcessConfigurator.h"
#include "FiniteDifference.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    ProcessConfigurator processConfigurator;
    processConfigurator.compute_process_split(256);
    std::cout << processConfigurator.y_proc << " " << processConfigurator.x_proc << std::endl;

    FiniteDifference(100, 100, 4, 3, 0, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 1, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 2, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 3, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 4, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 5, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 6, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 7, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 8, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 9, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 10, 0, 1, 0, 1);
    FiniteDifference(100, 100, 4, 3, 11, 0, 1, 0, 1);


    return 0;
}