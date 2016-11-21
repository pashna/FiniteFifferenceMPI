//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_FINITEDIFFERENCE_H
#define SUPERCOMPUTING_FINITEDIFFERENCE_H


class FiniteDifference {
private:
    // coordinates
    double x1, x2; // left and right border by x
    double y1, y2; // top and bottom border by y
    double hx, hy; // step x, step y;

    //
    int x_cell_start, y_cell_start;
    int x_cell_n, y_cell_n;

    bool left, right, top, bottom;

public: // actually private

    void compute_coordinates(int x_grid_n, int y_grid_n, int x_proc_n, int y_proc_n, int process_id);
};


#endif //SUPERCOMPUTING_FINITEDIFFERENCE_H
