//
// Created by pkochetk on 11/28/16.
//

#include "ResultWritter.h"

ResultWritter::ResultWritter(int grid_n, int proc_n, int proccess_id) {
    filename = "results/";
    std::string command = "mkdir -p ";
    filename += std::to_string(grid_n) + "x" + std::to_string(proc_n);
    command += filename;
    // create directory
    const int dir_err = system(command);

    if (-1 == dir_err)
    {
        std::cout << "Error creating directory!\n";
    }

    filename += std::to_string(proccess_id) + '.txt';
}

void ResultWritter::write_result(double *p, int x_cell_n, int y_cell_n) {
    std::ofstream result;
    result.open (filename);
    for (int j=0; j<y_cell_n; j++) {
        for (int i = 0; i < x_cell_n; i++)
            result << p[j * x_cell_n + i] << " ";
        result << "\n";
    }
    result.close();
}