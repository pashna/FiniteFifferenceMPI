//
// Created by pkochetk on 11/29/16.
//

#include "ProcessConfigurator.h"

void ProcessConfigurator::compute_process_split(int n_proc) {
    x_proc = 1;
    y_proc = 1;
    bool flag = true;
    while (x_proc * y_proc != n_proc) {
        if (flag)
            y_proc *= 2;
        else
            x_proc *= 2;

        flag = !flag;
    }
}