# FiniteDifference MPI + OpenMP

This algorithm solves the Dirichlet problem for Poisson's equation in a rectangular area using "steep descent iterations" for first several iterations and "conjugate gradient iterations" afterward.

####Parallel implementation of Dirichlet Problem for Poisson's Equation with Rectangular Area solved by Finite Difference Method.

Parallelism is provided by the bind of OpenMP and MPI technologies.
The program was tested on Lomonosov Supercomputer.

Method keypoints:

* five-point difference equation for Laplace operator approximation
* grid fragmentation are regular
* MPI technology for counting under supercomputers
* OpenMP technology for parallelization between processor cores scalar product: (a, b) = \sum_{i=1}^{i=n-1} ( \sum_{j=1}^{j=m-1} ( h'i * h'j * a(i, j) * b(i, j) )) 

Algorithm parameters: 
* boundary conditions: function 'fi' 
* right side of Laplace operator: function 'F' 
* stopping criteria: function 'stopCriteria'
