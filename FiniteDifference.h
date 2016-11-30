//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_FINITEDIFFERENCE_H
#define SUPERCOMPUTING_FINITEDIFFERENCE_H

#include "Condition.h"
#include "mpi.h"
#include <omp.h>
//#include <random>

class FiniteDifference {
public:
    FiniteDifference(Condition *_c, int x_grid_n, int y_grid_n, int x_proc_n_, int y_proc_n_, int process_id_,
                     double X1, double X2, double Y1, double Y2, MPI_Comm communicator_);

    int x_cell_n, y_cell_n;
    int iter;
    double *p;
    double norm;
    double hx, hy; // step x, step y;
    // coordinates for cell
    double x1, x2; // left and right border by x
    double y1, y2; // top and bottom border by y

private:


    // useful constants
    double hxhy, hx2, hy2;

    // cell coordinated
    int x_cell_start, y_cell_start;

    // border indicators
    bool left, right, top, bottom;

    int process_id, x_proc_n, y_proc_n;

    // algorithm
    double *p_prev, *delta_p;;
    double *g, *delta_g;
    double *r, *delta_r;

    Condition* c;

    // data to send and recieve
    double *send_lr, *receive_lr;
    double *send_rl, *receive_rl;
    double *send_td, *receive_td;
    double *send_bu, *receive_bu;
    MPI_Comm communicator;

public: // actually private

    // initialization phase
    void compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id);
    void initialize_constants(double X1, double X2, double Y1, double Y2, int x_grid_n, int y_grid_n);
    void initialize_matrixes();
    void initialize_mpi_arrays();

    // algorithm phase
    double max_norm();
    void compute_r();
    void compute_g(double alpha);
    void compute_p(double t);


    // tools
    double scalar_product(double* f1, double* f2);
    double drand(double fMin, double fMax) {
        //double f = (double)std::rand() / RAND_MAX;
        //return fMin + f * (fMax - fMin);
        return 0;
    }

    // finite difference computation
    void finite_difference(double* f, double* df);
    void write_MPI_arrays(double *f);
    void exchange_MPI_arrays();
    void difference_equation(double *f, double *df);
    void difference_equation_for_border(double *f, double *df);

public:
    void solve(double eps);
};


#endif //SUPERCOMPUTING_FINITEDIFFERENCE_H
