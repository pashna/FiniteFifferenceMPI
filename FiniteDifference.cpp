//
// Created by pkochetk on 11/20/16.
//

#include <iostream>
#include "FiniteDifference.h"
#include "Var6Cond.h"
#include <mpi.h>
#include "math.h"


FiniteDifference::FiniteDifference(Condition *_c, int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id,
                                   double X1, double X2, double Y1, double Y2) {
    compute_coordinates(x_grid_n, y_grid_n, x_proc_n, y_proc_n, process_id);
    initialize_constants(X1, X2, Y1, Y2, x_grid_n, y_grid_n);
    std::cout << "\n\n";
    c = _c;

}

void FiniteDifference::initialize_constants(double X1, double X2, double Y1, double Y2, int x_grid_n, int y_grid_n) {
    hx = (X2-X1)/(x_grid_n);
    hy = (Y2-Y1)/(y_grid_n);
    hxhy = hx*hy;
    hx2 = hx*hx;
    hy2 = hy*hy;

    x1 = x_cell_start * hx;
    x2 = (x_cell_start + x_cell_n) * hx;

    y1 = y_cell_start * hy;
    y2 = (y_cell_start + y_cell_n - 1) * hy;

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
        x_cell_start += (process_id % x_proc_n - (x_proc_n - x_proc_last_n)); // + offset

    }

    y_cell_start = process_id / x_proc_n * y_cell_n;

    std::cout << "!!! = " << y_proc_last_n << "\n";
    if (process_id / x_proc_n >= (y_proc_n - y_proc_last_n)) {
        // the process is charged for (Y_cells_n+1)  cells
        y_cell_n += 1;
        y_cell_start += (process_id / x_proc_n - (y_proc_n - y_proc_last_n)); // + offset
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


void FiniteDifference::initialize_matrixes() {
    int size = x_cell_n * y_cell_n;
    p = new double [size];
    p_prev = new double [size];
    delta_p = new double [size];
    g = new double [size];
    delta_g = new double [size];
    r = new double [size];
    delta_r = new double [size];

    double value;
    int index;
    if (top) {
        for(int i=0; i<x_cell_n; i++) {
            value = c->fi(x1 + i * hx, y1);
            index = i;
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + i * hx << ": " << y1 << "\t" << index << std::endl;

        }
    }
    if (right) {
        for(int i=0; i<y_cell_n; i++) {
            value = c->fi(x1 + (x_cell_n-1) * hx, y1 + hy * i);
            index = x_cell_n*(i+1)-1;
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + (x_cell_n-1) * hx << ": " << y1 + hy * i << "\t" << index << std::endl;

        }
    }
    if (bottom) {
        for(int i=0; i<x_cell_n; i++) {
            value = c->fi(x1 + i * hx, y2);
            index = x_cell_n*(y_cell_n-1) + i;
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + i * hx << ": " << y2 << "\t" << index << std::endl;
        }
    }
    if (left) {
        for(int i=0; i<y_cell_n; i++) {
            //std::cout << x1 << ": " << y1 + hy * i << std::endl;

            value = c->fi(x1, y1 + hy * i);
            index = x_cell_n*i;
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 << ": " << y1 + hy * i << "\t" << index << std::endl;
        }
    }


}


double FiniteDifference::scalar_product(double *f1, double *f2) {
    double s_local = 0;
    int index;
    // Start with y to avoid cache miss
    for (int j=0; j<y_cell_n; j++)
        for (int i=0; i<x_cell_n; i++) {
            index = j * x_cell_n + i;
            s_local += hxhy * f1[index] * f2[index];
        }

    double s_global = 0;

    // Summarize s from each of nodes and send to each of nodes
    MPI_Allreduce(&s_local, &s_global, 1, MPI_DOUBLE, MPI_SUM, procParams.comm);

    return s_global;

}


double FiniteDifference::max_norm() {
    double max_local = abs(p_prev[0]-p[0]);
    double tmp;
    int index;
    for (int j=0; j<y_cell_n; j++)
        for (int i=0; i<x_cell_n; i++) {
            index = j * x_cell_n + i;
            tmp = abs(p_prev[index] - p[index]);
            if (tmp > max_local)
                max_local = tmp;
        }

    double max_global;
    // Find max
    MPI_Allreduce(&max_local, &max_global, 1, MPI_DOUBLE, MPI_MAX, procParams.comm);
    return max_global;
}