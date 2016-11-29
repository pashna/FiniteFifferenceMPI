//
// Created by pkochetk on 11/28/16.
//

#ifndef SUPERCOMPUTING_RESULTWRITTER_H
#define SUPERCOMPUTING_RESULTWRITTER_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>



class ResultWritter {
public:
    ResultWritter(int grid_n, int proc_n, int proccess_id);
    void write_result(double* p, int x_cell_n, int y_cell_n, double x, double hx, double y, double hy);

private:
    std::string filename;
};


#endif //SUPERCOMPUTING_RESULTWRITTER_H
