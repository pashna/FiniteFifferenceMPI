#ifndef SUPERCOMPUTING_PROCESSCONFIGURATOR_H
#define SUPERCOMPUTING_PROCESSCONFIGURATOR_H

// The class decides how the grid will be shared between processes.
class ProcessConfigurator {
public:
    int x_proc;
    int y_proc;

    void compute_process_split(int n_proc);
};

#endif //SUPERCOMPUTING_PROCESSCONFIGURATOR_H
