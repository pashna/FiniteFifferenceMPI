//
// Created by pkochetk on 11/20/16.
//

#include <iostream>
#include "FiniteDifference.h"

void FiniteDifference::compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id) {
    // todo: easy to split for two different functions
    x_cell_n = (x_grid_n + 1) / x_proc_n;
    y_cell_n = (y_grid_n + 1) / y_proc_n;

    int x_cell_last_n = (x_grid_n + 1) - (x_proc_n * x_cell_n);
    int y_cell_last_n = (y_grid_n + 1) - (y_proc_n * y_cell_n);

    x_cell_start = (process_id % x_proc_n) * x_cell_n;
    if (process_id % x_proc_n >= (x_proc_n - x_cell_last_n)) {
        // the process is charged for (x_cells_n+1)  cells
        x_cell_n += 1;
        // + offset
        x_cell_start += (process_id % x_proc_n - (x_proc_n - x_cell_last_n));
    }

    y_cell_start = (process_id % x_proc_n) * y_cell_n;
    if (process_id % y_proc_n >= (y_proc_n - y_cell_last_n)) {
        // the process is charged for (x_cells_n+1)  cells
        y_cell_n += 1;
        // + offset
        y_cell_start += (process_id % y_proc_n - (y_proc_n - y_cell_last_n));
    }
    std::cout << process_id << std::endl;
    std::cout << "x_start="<<x_cell_start<<"; x_cell_n=" << x_cell_n << "\n";
    std::cout << "y_start="<<y_cell_start<<"; y_cell_n=" << y_cell_n << "\n";
    std::cout << "\n\n";

}
