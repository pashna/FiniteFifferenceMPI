#include <iostream>
#include "ProcessConfigurator.h"
#include "FiniteDifference.h"
int main() {
    std::cout << "Hello, World!" << std::endl;
    ProcessConfigurator processConfigurator;
    processConfigurator.compute_process_split(256);
    std::cout << processConfigurator.y_proc << " " << processConfigurator.x_proc << std::endl;

    FiniteDifference finiteDifference;
    finiteDifference.compute_coordinates(102, 102, 4, 4, 8);
    finiteDifference.compute_coordinates(102, 102, 4, 4, 9);
    finiteDifference.compute_coordinates(102, 102, 4, 4, 10);
    finiteDifference.compute_coordinates(102, 102, 4, 4, 11);

    return 0;
}