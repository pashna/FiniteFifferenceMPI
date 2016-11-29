#include <iostream>
#include "ProcessConfigurator.h"
#include "FiniteDifference.h"
#include "Var6Cond.h"
#include <mpi.h>
#include "ResultWritter.h"



int main(int argc, char** argv) {

    MPI_Init(&argc,&argv);
    Var6Cond cond;

    try {

        MPI_Comm communicator = MPI_COMM_WORLD;
        int process_id;
        int process_n;
        MPI_Comm_rank (communicator, &process_id); // get current process id
        MPI_Comm_size (communicator, &process_n); // get number of processes

        ProcessConfigurator proc_config;
        proc_config.compute_process_split(process_n);

        // Todo: move all soecific for the task params to cond and use cond as parameter for solve()
        FiniteDifference finiteDifference(&cond, 1000, 1000, proc_config.x_proc, proc_config.y_proc,
                                          process_id, 0, 1, 0, 1, communicator);
        finiteDifference.solve(0.0001);

        ResultWritter resultWritter(1000, process_n, process_id);
        resultWritter.write_result(finiteDifference.p, finiteDifference.x_cell_n, finiteDifference.y_cell_n);

    } catch (std::exception& e) {
        std::cout << e.what();
        std::cout << std::endl;
        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}