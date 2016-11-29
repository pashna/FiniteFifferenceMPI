#include <iostream>
#include "ProcessConfigurator.h"
#include "FiniteDifference.h"
#include "Var6Cond.h"
#include <mpi.h>
#include "ResultWritter.h"
#include <sys/time.h>



int main(int argc, char** argv) {

    int grid_n =  atoi(argv[1]);
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
        FiniteDifference finiteDifference(&cond, grid_n, grid_n, proc_config.x_proc, proc_config.y_proc,
                                          process_id, 0, 1, 0, 1, communicator);

        struct timeval tp;
        long int start_time = tp.tv_sec * 1000 + tp.tv_usec / 1000;

        if (process_id == 0) {
            gettimeofday(&tp, NULL);
        }

        finiteDifference.solve(0.0001);

        if (process_id == 0) {
            long int finish_time = tp.tv_sec * 1000 + tp.tv_usec / 1000;
            long int time = finish_time - start_time;
            std::cout << "=============\n";
            std::cout << "processes: " << process_n << "\n";
            std::cout << "grid: " << grid_n << "\n";
            std::cout << "time: " << time << "ms.\n";
            std::cout << "max_norm: " << finiteDifference.norm << "\n";
            std::cout << "iterations: " << finiteDifference.iter << "\n";
        }

        ResultWritter resultWritter(grid_n, process_n, process_id);
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