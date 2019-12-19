#include <vector>
#include <math.h>
#include <cassert>
#include <mpi.h>
#include <queue>
#include <algorithm>
#include "./MeshUtilities.h"

void parallelBucketSort(
	const std::vector<double3_t>& points_to_sort,
	std::vector<double3_t>& sorted_points,
	int sort_dim
);

void parallelMedian(
	const std::vector<double3_t>& local_points,
	int sort_dim,
	int P,
	double3_t& mdn_pt
);

bool sortDouble3tByX(const double3_t& p1, const double3_t& p2){ return p1.r[0] < p2.r[0];}
bool sortDouble3tByY(const double3_t& p1, const double3_t& p2){ return p1.r[1] < p2.r[1]; }
bool sortDouble3tByZ(const double3_t& p1, const double3_t& p2){ return p1.r[2] < p2.r[2]; }


//-----------------------------------------------------
// Usual Parallel Partitioning Code
//-----------------------------------------------------
void parallelRange(
	int globalstart,
	int globalstop,
	int irank,
	int nproc,
	int& localstart,
	int& localstop,
	int& localcount
)
{
	int nvals = globalstop - globalstart + 1;
	int divisor = nvals/nproc;
	int remainder = nvals%nproc;
	int offset;
	if (irank < remainder) offset = irank;
	else offset = remainder;

	localstart = irank*divisor + globalstart + offset;
	localstop = localstart + divisor - 1;
	if (remainder > irank) localstop += 1;
	localcount = localstop - localstart + 1;
}

//-----------------------------------------------------
// Extract a single coordinate list from a double3_t list
//-----------------------------------------------------
std::vector<double> getSingleCoordinateListFromPoints(
	const std::vector<double3_t>& points,
	int dim
)
{
    std::vector<double> coordinate_list(points.size());
    if (dim > 3 || dim < 0)
    {
        std::cerr << "Requested dimension " << dim << " is out of bounds at line " << __LINE__ << " of file " << __FILE__ << " for function " << __FUNCTION__ << std::endl;
        return coordinate_list;
    }

    for (unsigned int ipoint = 0; ipoint < points.size(); ipoint++)
    {
        coordinate_list[ipoint] = points[ipoint].r[dim];
    }

    //Want to see the list?
    /*
    for (unsigned int ipoint = 0; ipoint < coordinate_list.size(); ipoint++)
    {
        std::cout << coordinate_list[ipoint] << std::endl;
    }
    */
    return coordinate_list;
}


/*
This ORB function gets called by the other ORB function. This function
then recursively calls itself to perform ORB on a set of points,
assigning a label to each point.

This function works on a one-dimensional list of point labels, instead
of a two-dimensal list of local indices belonging to domains. Using
the label list makes the recursive point-labelling easier. Note that
while the number of unique labels will be P, some of the values will
be numbers greater than (P-1). The calling function is responsible
for sorting this out.
*/
void ORB(
	int P,
	const std::vector<double3_t>& points,
	std::vector<int>& labels,
	int dim
)
{

	int rank, nproc;
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(P>1)
	{
		int w0 = P/2;
		int w1 = P-w0;
		// Figure out what the global median is
		double3_t mdn_pt;
		parallelMedian(
			points,
			dim,
			P,
			mdn_pt
		);
		std::vector<double3_t> points_lt(points.size());
		std::vector<double3_t> points_geq(points.size());
		std::vector<int> lt_local_idcs(points.size());
		std::vector<int> geq_local_idcs(points.size());
		int n_lt_pts = 0,n_geq_pts = 0;
		// Bisect my points at the median.
		for(int ii = 0;ii<points.size();ii++)
		{
			if(points[ii].r[dim] < mdn_pt.r[dim])
			{
				lt_local_idcs[n_lt_pts] = ii;
				points_lt[n_lt_pts++] = points[ii];
			}
			else
			{
				geq_local_idcs[n_geq_pts] = ii;
				points_geq[n_geq_pts++] = points[ii];
			}
		}
		points_lt.resize(n_lt_pts);
		lt_local_idcs.resize(n_lt_pts);
		points_geq.resize(n_geq_pts);
		geq_local_idcs.resize(n_geq_pts);

		std::vector<int> lt_labels;
		std::vector<int> geq_labels;
		dim = (dim+1)%3;
		// ORB Strategy:
		// Bisect the supplied set of points into a left set and a
		// right set. Call this ORB function on each set. This will
		// create two lists of labels. They probably share some labels
		// i.e., the label '0' probably appears in lt_labels and in
		// geq_labels. In order to prevent these label collisions,
		// just map the left labels to even numbers (x2) and map the
		// right labels to odd numbers (x2+1).
		ORB(
			w0,
			points_lt,
			lt_labels,
			dim
		);
		ORB(
			w1,
			points_geq,
			geq_labels,
			dim
		);
		// Map the left labels to evens
		labels.resize(points.size());
		for(int ii = 0;ii<points_lt.size();ii++)
		{
			labels[lt_local_idcs[ii]] = 2*lt_labels[ii];
		}
		// Map the right labels to odds
		for(int ii = 0; ii < points_geq.size();ii++)
		{
			labels[geq_local_idcs[ii]] = 2*geq_labels[ii]+1;
		}
	}
	else
	{
		// In this implementation, every weight-1 set gets label 0.
		// Further up the calling tree, these labels will get changed.
		labels.resize(points.size());
		for(int ii = 0;ii<labels.size();ii++)
		{
			labels[ii]=0;
		}
	}
}


//-----------------------------------------------------
// Implementation of ORB.
//
// Inputs:
//
// P the number of domains to produce.
//
// points a list of type double3_t to partition.
//
// Output:
// points_in_orb_domains is a vector of vectors where
// points_in_orb_domains[iproc] stores the local element
// indexes of all processors in subdomain iproc (i.e.,
// destined for iproc).
// points_in_orb_domains[iproc][iele] = my local index of that point.
//-----------------------------------------------------
void ORB(
	int P,
	const std::vector<double3_t>& points,
	std::vector<std::vector<int> >& points_in_orb_domains
)
{
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	std::vector<int> labels;
	// Call the other ORB function (the one that does all the work).
	// Tell it to start bisection along the x dimension.
	ORB(
		nproc,
		points,
		labels,
		0
	);
	// The range of labels may be greater than nproc at this point.
	// Fix that.
    std::vector<int> pts_per_label(1<<(nproc-1));
	// Check which labels I have actually used.
    for(int ipt = 0; ipt < labels.size(); ipt++)
    {
        pts_per_label[labels[ipt]]++;
    }
    int label_to_check = 0;
    int cur_unique_label = 0;
    std::vector<int> label_map(1<<(nproc-1));
    while(cur_unique_label < nproc)
    {
        // Does any processor actually have points with this label?
        int n_pts_at_label = pts_per_label[label_to_check];
        MPI_Allreduce(
            &n_pts_at_label,
            &n_pts_at_label,
            1,
            MPI_INT,
            MPI_SUM,
            MPI_COMM_WORLD
        );

        // If there are points with this label, reassign them to
        // cur_unique_label
        if(n_pts_at_label)
		{
			label_map[label_to_check] = cur_unique_label++;
        }
        label_to_check++;
    }
	points_in_orb_domains.resize(nproc);
	for(int ivec = 0; ivec < points_in_orb_domains.size();ivec++)
	{
		points_in_orb_domains[ivec].resize(0);
	}
    for(int ipt = 0; ipt < labels.size(); ipt++)
    {
        labels[ipt] = label_map[labels[ipt]];
		points_in_orb_domains[labels[ipt]].push_back(ipt);
    }
}

void parallelMedian(
	const std::vector<double3_t>& local_points,
	int sort_dim,
	int P,
	double3_t& mdn_pt
)
{
	int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	std::vector<double3_t> sorted_points;

	// Perform the bucket sort.
	parallelBucketSort(
		local_points,
		sorted_points,
		sort_dim
	);

	// How many sorted points are on each proc?
	std::vector<int> sendcounts(nproc);
	std::vector<int> sdispls(nproc);
	std::vector<int> pts_per_proc(nproc);
	std::vector<int> recvcounts(nproc);
	std::vector<int> rdispls(nproc);
	int local_n_pts = sorted_points.size();
	for(int ii=0;ii<nproc;ii++)
	{
		sendcounts[ii] = 1;
		sdispls[ii] = 0;
		recvcounts[ii] = 1;
		rdispls[ii] = ii;
	}
	// Make number of points on each processor common knowledge.
	MPI_Alltoallv(
		&local_n_pts,
		&(sendcounts[0]),
		&(sdispls[0]),
		MPI_INT,
		&(pts_per_proc[0]),
		&(recvcounts[0]),
		&(rdispls[0]),
		MPI_INT,
		MPI_COMM_WORLD
	);
	int global_n_pts = 0;
	for(int ii = 0; ii < pts_per_proc.size();ii++)
	{
		global_n_pts += pts_per_proc[ii];
	}

	int w0 = P/2;
	int global_mdn_idx = ((global_n_pts)*w0)/P;
	int mdn_owner_local_idx = 0;
	int mdn_owner = 0;
	// Find out which processor has the median point.
	for(int ii = 0;ii<nproc;ii++)
	{
		if(global_mdn_idx >= pts_per_proc[ii])
		{
			global_mdn_idx -= pts_per_proc[ii];
		}
		else
		{
			mdn_owner = ii;
			mdn_owner_local_idx = global_mdn_idx;
			break;
		}
	}
	if(mdn_owner == rank)
	{
		mdn_pt = sorted_points[mdn_owner_local_idx];
	}
	// Broadcast the median.
	MPI_Bcast(
		&(mdn_pt),
		sizeof(double3_t),
		MPI_BYTE,
		mdn_owner,
		MPI_COMM_WORLD
	);
}


//-----------------------------------------------------
// Implementation of Parallel Bucket Sort.
//
// Inputs:
//
// The index of the dimension to sort along.
// The points to sort.
//
// Outputs:
//
// A subset of the sorted values. Note that the entries
// in sorted_points are not the same as points_to_sort
// on any given processor. They are a subset of the
// total set of sorted values where rank 0 will contain
// the lowest sorted numbers.
//-----------------------------------------------------
void parallelBucketSort(
	const std::vector<double3_t>& points_to_sort,
	std::vector<double3_t>& sorted_points,
	int sort_dim
)
{
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double minval = DBL_MAX;
	double maxval = -DBL_MAX;

	// Determine min and max values so that bucket sort can be
	// done with floor division. These need to be consistent across
	// all processors.
	for(int ii = 0;ii<points_to_sort.size();ii++)
	{

		if(points_to_sort[ii].r[sort_dim] < minval)
		{
			minval = points_to_sort[ii].r[sort_dim];
		}
		if(points_to_sort[ii].r[sort_dim] > maxval)
		{
			maxval = points_to_sort[ii].r[sort_dim];
		}
	}
	MPI_Allreduce(
		&minval,
		&minval,
		1,
		MPI_DOUBLE,
		MPI_MIN,
		MPI_COMM_WORLD
	);
	MPI_Allreduce(
		&maxval,
		&maxval,
		1,
		MPI_DOUBLE,
		MPI_MAX,
		MPI_COMM_WORLD
	);
	// Set up buckets
	std::vector<std::vector<double3_t> > send_buckets(nproc);
	// Place each point in the correct bucket.
	for(int ii = 0;ii<points_to_sort.size();ii++)
	{
		double ival = points_to_sort[ii].r[sort_dim];
		int b_idx = (int)((ival-minval)/(maxval-minval)*nproc);
		if(b_idx >= nproc) b_idx = nproc - 1;
		send_buckets[b_idx].push_back(points_to_sort[ii]);
	}

	std::vector<std::vector<double3_t> > recv_buckets(nproc);
	// Swap Buckets
	MPI_Alltoall_vecvecT(send_buckets,recv_buckets);
	// Now I have a bunch of buckets from different processors.
	// Put them all in one bucket.
	// First, determine how big that bucket needs to be
	int n_recv_pts = 0;
	for(int bb = 0;bb<recv_buckets.size();bb++)
	{
		n_recv_pts += recv_buckets[bb].size();
	}
	sorted_points.resize(n_recv_pts);
	// Now dump all the small buckets into the big bucket.
	int sorted_pt_idx = 0;
	for(int bb = 0;bb<recv_buckets.size();bb++)
	{
		for(int ii = 0; ii < recv_buckets[bb].size();ii++)
		{
			sorted_points[sorted_pt_idx++] = recv_buckets[bb][ii];
		}
	}
	// Now, sort the big bucket.
	switch(sort_dim)
	{
		default:
		case 0:
			sort(sorted_points.begin(), sorted_points.end(),sortDouble3tByX);
			break;
		case 1:
			sort(sorted_points.begin(), sorted_points.end(),sortDouble3tByY);
			break;
		case 2:
			sort(sorted_points.begin(), sorted_points.end(),sortDouble3tByZ);
			break;
	}
}

void qsortPoints(
	std::vector<double3_t>& points,
	int dim
)
{
	if(points.size() < 2) return;

	double3_t pivot = points[0];
	//std::cout << "Pivot: " << pivot.r[dim] << std::endl;
	std::vector<double3_t> lesser,greater;
	lesser.resize(points.size());
	greater.resize(points.size());
	int l_idx=0,g_idx=0;
	for(int ii = 1;ii<points.size();ii++)
	{
		if(points[ii].r[dim] < pivot.r[dim])
		{
			lesser[l_idx++] = points[ii];
			//std::cout << "Lesser gets " << points[ii].r[dim] << std::endl;
		}
		else
		{
			greater[g_idx++]=points[ii];
			//std::cout << "Greater gets " << points[ii].r[dim] << std::endl;
		}
	}
	lesser.resize(l_idx);
	greater.resize(g_idx);

	qsortPoints(lesser,dim);
	qsortPoints(greater,dim);
	int p_idx = 0;
	for(int ii = 0; ii < lesser.size(); ii++)
	{
		points[p_idx++] = lesser[ii];
	}
	points[p_idx++] = pivot;
	for(int ii = 0; ii < greater.size(); ii++)
	{
		points[p_idx++] = greater[ii];
	}
}
