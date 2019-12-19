#include <iostream>
#include <mpi.h>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include "Mesh.h"




int main(int argc, char** argv)
{
    int rank,nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    srand(MPI_Wtime() + rank);

    std::stringstream mesh_name;
    std::stringstream cell_filename;
    cell_filename << "CellMatrix_P" << nproc << ".txt";
    mesh_name << argv[1] << ".msh";

    int n_iter = atoi(argv[2]);
    double pop_percent = atof(argv[3]);
    double pop_rate = (pop_percent)/100.0;
    bool log_each_iteration = (atoi(argv[4])>0);
    std::cout << "pop rate = " << pop_rate << std::endl;
    Mesh BallMesh(mesh_name.str());

    BallMesh.partitionMesh();
    BallMesh.calculateVertexConnectivity();
    BallMesh.populateMeshVertices(pop_rate);
    BallMesh.writeCellMatrixFile(cell_filename.str(),false);
    //BallMesh.outputStatistics();
    BallMesh.updateVertexStates();
    double tstart = MPI_Wtime();
    double tnow;
    double iter_per_s;
    MPI_Barrier(MPI_COMM_WORLD);
    int disp_counter = 0;
    for(int iter = 0; iter < n_iter; iter++)
    {
        if(log_each_iteration)
        {
            BallMesh.writeCellMatrixFile(cell_filename.str(),true);
        }
        BallMesh.updateVertexStates();
        if((disp_counter++ > 200))
        {
            if(rank == 0)
            {
                disp_counter = 0;
                tnow = MPI_Wtime();
                iter_per_s = iter/(tnow-tstart);
                std::cout << "\r" << iter << " / " << n_iter << ", " << iter_per_s << std::flush;
            }
        }
    }
    if(rank == 0) std::cout << std::endl;
    BallMesh.writeCellMatrixFile(cell_filename.str(),true);
    MPI_Finalize();

    return(0);
}
