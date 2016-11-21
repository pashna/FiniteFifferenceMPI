//
// Created by pkochetk on 11/20/16.
//

#ifndef SUPERCOMPUTING_PROCESSCONFIGURATOR_H
#define SUPERCOMPUTING_PROCESSCONFIGURATOR_H


class ProcessConfigurator {
    // The class decides how the grid will be shared between processes.
public:
    int x_proc;
    int y_proc;

    void compute_process_split(int n_proc) {
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

    };
};


#endif //SUPERCOMPUTING_PROCESSCONFIGURATOR_H
