//
// Created by pkochetk on 11/20/16.
//

#include <iostream>
#include "FiniteDifference.h"

FiniteDifference::FiniteDifference(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id,
                                   double X1, double X2, double Y1, double Y2) {
    compute_coordinates(x_grid_n, y_grid_n, x_proc_n, y_proc_n, process_id);
    initialize_constants(X1, X2, Y1, Y2, x_grid_n, y_grid_n);
    std::cout << "\n\n";

}

void FiniteDifference::initialize_constants(double X1, double X2, double Y1, double Y2, int x_grid_n, int y_grid_n) {
    hx = (X2-X1)/(x_grid_n+1);
    hy = (Y2-Y1)/(y_grid_n+1);
    hxhy = hx*hy;
    hx2 = hx*hx;
    hy2 = hy*hy;

    x1 = x_cell_start * hx;
    x2 = (x_cell_start + x_cell_n) * hx;

    y1 = y_cell_start * hy;
    y2 = (y_cell_start + y_cell_n) * hy;

    std::cout << "y:(" << y1 << "; " << y2 << ")  x:(" << x1 << "; " << x2 << ")";

}


void FiniteDifference::compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id) {
    // todo: easy to split for two different functions
    x_cell_n = (x_grid_n + 1) / x_proc_n;
    y_cell_n = (y_grid_n + 1) / y_proc_n;

    int x_proc_last_n = (x_grid_n + 1) - (x_proc_n * x_cell_n);
    int y_proc_last_n = (y_grid_n + 1) - (y_proc_n * y_cell_n);

    x_cell_start = (process_id % x_proc_n) * x_cell_n;
    if (process_id % x_proc_n >= (x_proc_n - x_proc_last_n)) {
        // the process is charged for (x_cells_n+1)  cells
        x_cell_n += 1;
        // + offset
        x_cell_start += (process_id % x_proc_n - (x_proc_n - x_proc_last_n));
    }

    y_cell_start = process_id / x_proc_n * y_cell_n;

    std::cout << "!!! = " << y_proc_last_n << "\n";
    if (process_id / x_proc_n >= (y_proc_n - y_proc_last_n)) {
        // the process is charged for (Y_cells_n+1)  cells
        y_cell_n += 1;
        // + offset
        y_cell_start += (process_id / x_proc_n - (y_proc_n - y_proc_last_n));
    }

    if (process_id % x_proc_n == 0)
        left = true;
    else
        left = false;

    if (process_id % x_proc_n == x_proc_n -1)
        right = true;
    else
        right = false;

    if (process_id <= x_proc_n)
        top = true;
    else
        top = false;

    if (process_id >= x_proc_n * (y_proc_n - 1))
        bottom = true;
    else
        bottom = false;

    std::cout << process_id << std::endl;
    std::cout << "x_start="<<x_cell_start<<"; x_cell_n=" << x_cell_n << "\n";
    std::cout << "y_start="<<y_cell_start<<"; y_cell_n=" << y_cell_n << "\n";
    std::cout << "t="<<top << " r=" << right << " b=" << bottom << " l=" << left << "\n";
}
