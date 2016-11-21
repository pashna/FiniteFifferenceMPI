//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_FINITEDIFFERENCE_H
#define SUPERCOMPUTING_FINITEDIFFERENCE_H

#include "Condition.h"


class FiniteDifference {
public:
    FiniteDifference(Condition *_c, int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id,
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

    // algorithm
    double* p;
    double* p_prev;

    Condition* c;

public: // actually private

    void compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id);
    void initialize_constants(double X1, double X2, double Y1, double Y2, int x_grid_n, int y_grid_n);
    double scalar_product(double* f1, double* f2);
    void initialize_matrixes();
};


#endif //SUPERCOMPUTING_FINITEDIFFERENCE_H
