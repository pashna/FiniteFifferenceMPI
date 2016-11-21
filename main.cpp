#include <iostream>
#include "ProcessConfigurator.h"
#include "FiniteDifference.h"
#include "Var6Cond.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    ProcessConfigurator processConfigurator;
    processConfigurator.compute_process_split(256);
    std::cout << processConfigurator.y_proc << " " << processConfigurator.x_proc << std::endl;

    Var6Cond cond;
    FiniteDifference a(&cond, 100, 100, 4, 3, 0, 0, 1, 0, 1);
    FiniteDifference a2(&cond, 100, 100, 4, 3, 4, 0, 1, 0, 1);
    FiniteDifference a3(&cond, 100, 100, 4, 3, 8, 0, 1, 0, 1);


    a.initialize_matrixes();
    std::cout << "\n";
    a2.initialize_matrixes();
    std::cout << "\n";
    a3.initialize_matrixes();
    std::cout << "\n";



    return 0;
}