//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_FINITEDIFFERENCE_H
#define SUPERCOMPUTING_FINITEDIFFERENCE_H

#include "Condition.h"


class FiniteDifference {
public:
    FiniteDifference(Condition *_c, int x_grid_n, int y_grid_n, int x_proc_n_, int y_proc_n_, int process_id_,
                     double X1, double X2, double Y1, double Y2);
private:
    // coordinates for cell
    double x1, x2; // left and right border by x
    double y1, y2; // top and bottom border by y
    double hx, hy; // step x, step y;

    // useful constants
    double hxhy, hx2, hy2;

    // cell coordinated
    int x_cell_start, y_cell_start;
    int x_cell_n, y_cell_n;

    // border indicators
    bool left, right, top, bottom;

    int process_id, x_proc_n, y_proc_n;

    // algorithm
    double* p;
    double* p_prev;
    double* g;
    double* r;
    double* delta_p;
    double* delta_g;
    double* delta_r;

    Condition* c;

    // data to send and recieve
    double* send_lr;
    double* receive_lr;
    double* send_rl;
    double* receive_rl;
    double* send_td;
    double* receive_td;
    double* send_bu;
    double* receive_bu;

public: // actually private

    void compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id);
    void initialize_constants(double X1, double X2, double Y1, double Y2, int x_grid_n, int y_grid_n);
    void initialize_matrixes();
    void initialize_mpi_arrays();
    double scalar_product(double* f1, double* f2);
    double max_norm();
    double drand(double fMin, double fMax) {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }
    void solve(double eps);

    // finite difference
    void finite_difference(double* f, double* df);
    void write_MPI_arrays(double *f);
    void exchange_MPI_arrays();
    void difference_equation(double *f, double *df);
    void difference_equation_for_border(double *f, double *df);

};


#endif //SUPERCOMPUTING_FINITEDIFFERENCE_H
