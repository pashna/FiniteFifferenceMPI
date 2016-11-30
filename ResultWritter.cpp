//
// Created by pkochetk on 11/28/16.
//

#include "ResultWritter.h"

ResultWritter::ResultWritter(int grid_n, int proc_n, int proccess_id) {
    filename = "results/";
    std::string command = "mkdir -p ";
    filename += number_to_string(grid_n) + "x" + number_to_string(proc_n);
    command += filename;
    // create directory
    int dir_err = system(command.c_str());

    if (-1 == dir_err) {
        std::cout << "Error creating directory!\n";
    }
    filename += number_to_string(proccess_id) + ".txt";
}

void ResultWritter::write_result(double *p, int x_cell_n, int y_cell_n, double x, double hx, double y, double hy) {
    std::ofstream result;
    result.open (filename.c_str());
    for (int j=0; j<y_cell_n; j++) {
        for (int i = 0; i < x_cell_n; i++) {
            result << p[j * x_cell_n + i] << " " << x + hx * i << " " << y + hy * j << "\t";

        }
        result << "\n";
    }
    result.close();
}
