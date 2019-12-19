#ifndef _MESH_UTILITIES_H_
#define _MESH_UTILITIES_H_

#include <vector>
#include <iostream>
#include <cfloat>

typedef struct
{
    double r[3];
}
double3_t;

std::vector<double> getSingleCoordinateListFromPoints(const std::vector<double3_t>& points, int dim);

//You are responsible for implementing both of these functions in MeshUtilities.cpp.
void ORB(
    int P,
    const std::vector<double3_t>& points,
    std::vector<std::vector<int> >& points_in_orb_domains
);
void ORB(
	int P,
	const std::vector<double3_t>& points,
	std::vector<int>& labels,
	int dim
);
void parallelBucketSort(const std::vector<double>& values_to_sort, std::vector<double>& sorted_values);
void qsortPoints(std::vector<double3_t>& points,int dim);
void parallelMedian(
	const std::vector<double3_t>& local_points,
	int sort_dim,
	int P,
	double3_t& mdn_pt
);


//Usual parallel range function
void parallelRange(int globalstart, int globalstop, int irank, int nproc, int& localstart, int& localstop, int& localcount);

//--------------------------------------------------------------------------
// Utility function for passing a vector of vectors of some type T in an
// all-to-all variable fashion. Data to recv is a vector of vector of type T
// where data_to_recv[iproc] is data received from iproc.
//--------------------------------------------------------------------------
template<class T>
void MPI_Alltoall_vecvecT(const std::vector<std::vector<T> >& data_to_send, std::vector<std::vector<T> >& data_to_recv)
{
    int nproc;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    assert(nproc == (int)data_to_send.size());

    std::vector<int> send_counts(nproc);
    std::vector<int> recv_counts(nproc);
    std::vector<int> send_displacements(nproc);
    std::vector<int> recv_displacements(nproc);

    int total_send_count = 0; //in bytes
    
    for (unsigned int iproc = 0; iproc < nproc; iproc++)
    {
        send_counts[iproc] = data_to_send[iproc].size()*sizeof(T);
        send_displacements[iproc] = total_send_count;
        total_send_count += send_counts[iproc];
    }

    MPI_Alltoall(&send_counts[0],1,MPI_INT, &recv_counts[0], 1, MPI_INT, MPI_COMM_WORLD);

    int total_recv_count = 0;   //in bytes
    for (unsigned int iproc = 0; iproc < nproc; iproc++)
    {
        recv_displacements[iproc] = total_recv_count;
        total_recv_count += recv_counts[iproc];
    }

    assert(total_send_count%sizeof(T)==0);
    assert(total_recv_count%sizeof(T)==0);

    std::vector<T> sendbuf(total_send_count/sizeof(T));
    std::vector<T> recvbuf(total_recv_count/sizeof(T));


    int count = 0;
    for (int i = 0; i < nproc; i++)
    {
        for (unsigned j = 0; j < data_to_send[i].size(); j++)
        {
            sendbuf[count] = data_to_send[i][j];
            count++;
        }
    }

    MPI_Alltoallv(&sendbuf[0],&send_counts[0],&send_displacements[0],MPI_BYTE,
				  &recvbuf[0],&recv_counts[0],&recv_displacements[0],MPI_BYTE,
				  MPI_COMM_WORLD);

    data_to_recv.resize(nproc);
    count = 0;
    for (int i = 0; i < nproc; i++)
    {
        data_to_recv[i].resize(0);
        for (int j = 0; j < recv_counts[i]/sizeof(T); j++)  //remember recv counts are in bytes
        {
        	data_to_recv[i].push_back(recvbuf[count]);
            count++;
        }
    }

}

#endif
