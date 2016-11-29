//
// Created by pkochetk on 11/20/16.
//

#include <iostream>
#include "FiniteDifference.h"
#include "Var6Cond.h"
#include "math.h"
//#include <random>

FiniteDifference::FiniteDifference(Condition *_c, int x_grid_n, int y_grid_n, int x_proc_n_, int y_proc_n_, int process_id_,
                                   double X1, double X2, double Y1, double Y2, MPI_Comm communicator_) {
    compute_coordinates(x_grid_n, y_grid_n, x_proc_n_, y_proc_n_, process_id_);
    initialize_constants(X1, X2, Y1, Y2, x_grid_n, y_grid_n);
    initialize_mpi_arrays();
    c = _c;
    process_id = process_id_;
    communicator = communicator_;
}

void FiniteDifference::initialize_mpi_arrays() {
        send_lr = new double [y_cell_n];
        receive_lr = new double [y_cell_n];

        send_rl = new double [y_cell_n];
        receive_rl = new double [y_cell_n];

        send_td = new double [x_cell_n];
        receive_td = new double [x_cell_n];

        send_bu = new double [x_cell_n];
        receive_bu = new double [x_cell_n];
}

void FiniteDifference::solve(double eps) {
    initialize_matrixes();
    double scalar_dg_g = 1;
    double scalar_dr_g = 1;
    double scalar_r_g = 1;
    double a = 0;
    double t = 0;

    // first iteration
    finite_difference(delta_p, p_prev);
    compute_r();
    std::swap(g, r);
    finite_difference(delta_g, g);
    scalar_r_g = scalar_product(g, g);
    scalar_dg_g = scalar_product(delta_g, g);
    t = scalar_r_g / scalar_dg_g;
    compute_p (t);
    norm = max_norm();
    while(norm >= eps) {
        std::swap(p, p_prev);
        // each row is a step of the algorithm
        finite_difference(delta_p, p_prev);
        compute_r();
        finite_difference(r, delta_p);
        finite_difference(delta_r, r);
        scalar_dr_g = scalar_product(delta_r, g);
        a = scalar_dr_g / scalar_dg_g;
        compute_g(a);
        finite_difference(delta_g, g);
        scalar_r_g = scalar_product(r, g);
        scalar_dg_g = scalar_product(delta_g, g);
        t = scalar_r_g / scalar_dg_g;
        compute_p (t);
        norm = max_norm();
        iter++;
    } ;
}

void FiniteDifference::finite_difference(double *f, double *df) {
    // calculate difference equation for each inner node
    difference_equation(f, df);
    write_MPI_arrays(f);
    exchange_MPI_arrays();
    difference_equation_for_border(f, df);
}

void FiniteDifference::difference_equation_for_border(double *f, double *df) {
    int cur_idx;
    int left_idx;
    int right_idx;
    int up_idx;
    int down_idx;

    // left->right
    int i=0;
    int j=0;
    if (!left) {
        for (j = 1; j < y_cell_n - 1; j++) {
            cur_idx = j * x_cell_n + i;
            left_idx = cur_idx-1;
            right_idx = cur_idx+1;
            up_idx = (j - 1) * x_cell_n + i;
            down_idx = (j + 1) * x_cell_n + i;

            df[cur_idx] = (2*f[cur_idx] - receive_lr[j] - f[right_idx])/hx2 +
                       (2*f[cur_idx] - f[up_idx] - f[down_idx])/hy2;
        }
    }
    if (!top & !left) {
        df[0] = (2*f[0] - receive_lr[0] - f[1])/hx2 +
                (2*f[0] - receive_td[0] - f[x_cell_n])/hy2;
    }


    // right->left
    i=x_cell_n-1;;
    j=0;
    if (!right) {
        for (j = 1; j < y_cell_n - 1; j++) {
            cur_idx = j * x_cell_n + i;
            left_idx = cur_idx-1;
            right_idx = cur_idx+1;
            up_idx = (j - 1) * x_cell_n + i;
            down_idx = (j + 1) * x_cell_n + i;

            df[cur_idx] = (2*f[cur_idx] - f[left_idx]- receive_rl[j])/hx2 +
                          (2*f[cur_idx] - f[up_idx] - f[down_idx])/hy2;
        }
    }
    if (!top & !right) {
        df[x_cell_n-1] = (2*f[x_cell_n-1] - f[x_cell_n-2] - receive_rl[0] )/hx2 +
                         (2*f[x_cell_n-1] - receive_td[x_cell_n-1] - f[2*x_cell_n-1])/hy2;
    }


    // top -> down
    i=0;
    j=0;
    if (!top) {
        for (i = 1; i < x_cell_n - 1; i++) {
            cur_idx = j * x_cell_n + i;
            left_idx = cur_idx-1;
            right_idx = cur_idx+1;
            up_idx = (j - 1) * x_cell_n + i;
            down_idx = (j + 1) * x_cell_n + i;

            df[cur_idx] = (2*f[cur_idx] - f[left_idx] - f[right_idx])/hx2 +
                          (2*f[cur_idx] - f[up_idx] - receive_td[i])/hy2;
        }
    }
    if (!bottom & !left){
        df[(y_cell_n) * (x_cell_n-1)] = (2*f[(y_cell_n) * (x_cell_n-1)] - receive_lr[y_cell_n -1] - receive_rl[(y_cell_n) * (x_cell_n-1) + 1] )/hx2 +
                                        (2*f[(y_cell_n) * (x_cell_n-1)] - f[(y_cell_n) * (x_cell_n-2)] - receive_bu[0])/hy2;
    }


    // down -> top
    i=0;
    j=y_cell_n -1;
    if (!bottom) {
        for (i = 1; i < x_cell_n - 1; i++) {
            cur_idx = j * x_cell_n + i;
            left_idx = cur_idx-1;
            right_idx = cur_idx+1;
            up_idx = (j - 1) * x_cell_n + i;
            down_idx = (j + 1) * x_cell_n + i;

            df[cur_idx] = (2*f[cur_idx] - f[left_idx] - f[right_idx])/hx2 +
                          (2*f[cur_idx] - receive_bu[i] - f[down_idx])/hy2;
        }
    }
    if (!bottom & !right){
        cur_idx = y_cell_n * x_cell_n - 1 ;
        df[cur_idx] = (2*f[cur_idx] - receive_lr[cur_idx - 1] - receive_rl[y_cell_n -1] )/hx2 +
                      (2*f[cur_idx] - f[cur_idx - x_cell_n] - receive_bu[0])/hy2;
    }


}

void FiniteDifference::exchange_MPI_arrays() {
    /*
     * TODO: TAGS!
     */
    //left -> right
    int ret;

    if (!right) {
        ret = MPI_Send(send_lr, y_cell_n, MPI_DOUBLE, process_id + 1, 0, communicator);
        if (ret != MPI_SUCCESS) throw 5;
    }
    if (!left) {
        ret = MPI_Recv(receive_lr, y_cell_n, MPI_DOUBLE, process_id - 1, 0, communicator, MPI_STATUS_IGNORE);
        if (ret != MPI_SUCCESS) throw 5;
    }
    // right -> left
    if (!left) {
        ret = MPI_Send(send_rl, y_cell_n, MPI_DOUBLE, process_id - 1, 0, communicator);
        if (ret != MPI_SUCCESS) throw 5;
    }
    if (!right) {
       ret = MPI_Recv(receive_rl, y_cell_n, MPI_DOUBLE, process_id + 1, 0, communicator, MPI_STATUS_IGNORE);
        if (ret != MPI_SUCCESS) throw 5;
    }
    // bottom -> up
    if (!top) {
        ret = MPI_Send(send_td, x_cell_n, MPI_DOUBLE, process_id + x_proc_n, 0, communicator);
        if (ret != MPI_SUCCESS) throw 5;
    }
    if (!bottom) {
        ret = MPI_Recv(receive_td, x_cell_n, MPI_DOUBLE, process_id - x_proc_n, 0, communicator, MPI_STATUS_IGNORE);
        if (ret != MPI_SUCCESS) throw 5;
    }

    // up -> bottom
    if (!top) {
        ret = MPI_Send(send_bu, x_cell_n, MPI_DOUBLE, process_id - x_proc_n, 0, communicator);
        if (ret != MPI_SUCCESS) throw 5;
    }
    if (!bottom) {
        ret = MPI_Recv(receive_bu, x_cell_n, MPI_DOUBLE, process_id + x_proc_n, 0, communicator, MPI_STATUS_IGNORE);
        if (ret != MPI_SUCCESS) throw 5;
    }

}

void FiniteDifference::difference_equation(double *f, double *df) {
    int cur_idx;
    int left_idx;
    int right_idx;
    int up_idx;
    int down_idx;

    #pragma omp parallel
    #pragma omp for schedule (static)
    for (int j = 1; j < y_cell_n - 1; j++)
        for (int i = 1; i < x_cell_n - 1; i++) {
            cur_idx = j * x_cell_n + i;
            left_idx = cur_idx - 1;
            right_idx = cur_idx + 1;
            up_idx = (j - 1) * x_cell_n + i;
            down_idx = (j + 1) * x_cell_n + i;

            df[cur_idx] = (2 * f[cur_idx] - f[left_idx] - f[right_idx]) / hx2 +
                          (2 * f[cur_idx] - f[up_idx] - f[down_idx]) / hy2;
        }

}

void FiniteDifference::write_MPI_arrays(double *f) {
    // bottom -> up
    int index;
    for(int i=0; i<x_cell_n; i++) {
        index = i;
        send_bu[i] = f[index];
    }

    // left->right
    for(int i=0; i<y_cell_n; i++) {
        index = x_cell_n * (i + 1) - 1;
        send_lr[i] = f[index];
    }

    // top->bottom
    for(int i=0; i<x_cell_n; i++) {
        index = x_cell_n * (y_cell_n - 1) + i;
        send_td[i] = f[index];
    }

    // right->left
    for(int i=0; i<y_cell_n; i++) {
        index = x_cell_n * i;
        send_rl[i] = f[index];
    }


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

    //std::cout << "y:(" << y1 << "; " << y2 << ")  x:(" << x1 << "; " << x2 << ")";

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

    //std::cout << process_id << std::endl;
    //std::cout << "x_start="<<x_cell_start<<"; x_cell_n=" << x_cell_n << "\n";
    //std::cout << "y_start="<<y_cell_start<<"; y_cell_n=" << y_cell_n << "\n";
    //std::cout << "t="<<top << " r=" << right << " b=" << bottom << " l=" << left << "\n";
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

    for (int j=0; j<y_cell_n; j++)
        for (int i=0; i<x_cell_n; i++)
            p_prev[j * x_cell_n + i] = drand(-100, 100);


    if (top) {
        for(int i=0; i<x_cell_n; i++) {
            index = i;
            value = c->fi(x1 + i * hx, y1);
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + i * hx << ": " << y1 << "\t" << index << std::endl;

        }
    }
    if (right) {
        for(int i=0; i<y_cell_n; i++) {
            index = x_cell_n*(i+1)-1;
            value = c->fi(x1 + (x_cell_n-1) * hx, y1 + hy * i);
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + (x_cell_n-1) * hx << ": " << y1 + hy * i << "\t" << index << std::endl;

        }
    }
    if (bottom) {
        for(int i=0; i<x_cell_n; i++) {
            index = x_cell_n*(y_cell_n-1) + i;
            value = c->fi(x1 + i * hx, y2);
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 + i * hx << ": " << y2 << "\t" << index << std::endl;
        }
    }
    if (left) {
        for(int i=0; i<y_cell_n; i++) {
            index = x_cell_n*i;
            value = c->fi(x1, y1 + hy * i);
            p_prev[index] = value;
            p[index] = value;
            delta_p[index] = g[index] = delta_g[index] = r[index] = delta_r[index] = 0;
            //std::cout << x1 << ": " << y1 + hy * i << "\t" << index << std::endl;
        }
    }


}

void FiniteDifference::compute_r() {
    int index;
    // biases
    int t_ = (top) ? 1 : 0;
    int b_ = (bottom) ? 1 : 0;
    int r_ = (right) ? 1 : 0;
    int l_ = (left) ? 1 : 0;

    for (int j = t_; j < y_cell_n - b_; j++){
        for (int i = l_; i < x_cell_n - r_; i++){
            index = j * x_cell_n + i;
            r[index] = delta_p[index] - c->F(x1 + i * hx, y2 + j * hy);
        }
    }
}

void FiniteDifference::compute_g(double a) {
    int index;
    // biases
    int t_ = (top) ? 1 : 0;
    int b_ = (bottom) ? 1 : 0;
    int r_ = (right) ? 1 : 0;
    int l_ = (left) ? 1 : 0;

    for (int j = t_; j < y_cell_n - b_; j++){
        for (int i = l_; i < x_cell_n - r_; i++){
            index = j * x_cell_n + i;
            g[index] = r[index] - a * g[index];
        }
    }
}

void FiniteDifference::compute_p(double t) {
    int index;
    // biases
    int t_ = (top) ? 1 : 0;
    int b_ = (bottom) ? 1 : 0;
    int r_ = (right) ? 1 : 0;
    int l_ = (left) ? 1 : 0;

    for (int j = t_; j < y_cell_n - b_; j++){
        for (int i = l_; i < x_cell_n - r_; i++){
            index = j * x_cell_n + i;
            p[index] = p_prev[index] - t * g[index];
        }
    }
}

double FiniteDifference::scalar_product(double *f1, double *f2) {
    double s_local = 0;
    int index;
    #pragma omp parallel
    #pragma omp for schedule(static) reduction(+:s_local)
    // Start with y to avoid cache miss
    for (int j = 0; j < y_cell_n; j++)
        for (int i = 0; i < x_cell_n; i++) {
            index = j * x_cell_n + i;
            s_local += hxhy * f1[index] * f2[index];
        }

    double s_global = 0;

    // Summarize s from each of nodes and send to each of nodes
    MPI_Allreduce(&s_local, &s_global, 1, MPI_DOUBLE, MPI_SUM, communicator);

    return s_global;

}


double FiniteDifference::max_norm() {
    double max_local = fabs(p_prev[0]-p[0]);
    double tmp;
    int index;
    for (int j=0; j<y_cell_n; j++)
        for (int i=0; i<x_cell_n; i++) {
            index = j * x_cell_n + i;
            tmp = fabs(p_prev[index] - p[index]);
            if (tmp > max_local)
                max_local = tmp;
        }

    double max_global;
    // Find max
    MPI_Allreduce(&max_local, &max_global, 1, MPI_DOUBLE, MPI_MAX, communicator);
    return max_global;
}