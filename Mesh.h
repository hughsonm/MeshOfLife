#ifndef _MESH_H_
#define _MESH_H_

#include <iostream>
#include <vector>
#include <limits.h>
#include <fstream>
#include <math.h>
#include <cassert>
#include <mpi.h>
#include <map>

#include <queue>
#include <sstream>
#include "./MeshUtilities.h"

#define MAX_ELEMENT_VERTS 4

enum FV_MESH_GMESH_ELEMENT
{
	FV_MESH_MIN_GMESH_ELEMENT,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_LINE = 1,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_TRIANGLE = 2,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_QUADRANGLE = 3,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_TETRAHEDRAL = 4,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_HEXAHEDRAL = 5,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_PRISM = 6,
	FV_MESH_GMESH_ELEMENT_FIRST_ORDER_PYRAMID = 7,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_LINE = 8,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_TRIANGLE = 9,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_QUADRANGLE = 10,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_TETRAHEDRAL = 11,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_HEXAHEDRAL = 12,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_PRISM = 13,
	FV_MESH_GMESH_ELEMENT_SECOND_ORDER_PYRAMID = 14,
	FV_MESH_GMESH_ELEMENT_POINT = 15,
	FV_MESH_MAX_GMESH_ELEMENT
};


//using namespace std;

typedef struct helper_vertex_t
{
	int host_mid;
	int nbr_mid;
	bool nbr_state;
}
helper_vertex_t;

typedef struct vertex_t
{
    double r[3];						// Position in 3-space
    int mid;							// Mesh ID
	bool current_state;					// Alive or dead
	bool next_state;					// Alive or dead
	int owner;							// Rank who owns vertex
	bool is_shared;						// Tell me whether this is shared
	int alive_dead[2];					// Number of alive/dead neighbours
	std::vector<int> family;			// List of ranks which share vertex
	std::vector<vertex_t*> neighbours;	// List of pointers to local neighbouring vertices
	std::vector<int> global_nbr_mids;	// List of global connections
	std::vector<bool> nbr_was_counted;	// Used during vertex update
}
vertex_t;

typedef struct element_t
{
    int vertex_mids[MAX_ELEMENT_VERTS];
    int nvert;
    int mid;
    int type;
	bool current_state;
	bool next_state;
	element_t * neighbours[4];
}
element_t;


inline bool sortVertexByMIDPredicate(const vertex_t& v1, const vertex_t& v2){ return v1.mid < v2.mid; }
inline int searchVertexByMIDPredicate(const void* v1, const void* v2) {return ((const vertex_t*)v1)->mid - ((const vertex_t*)v2)->mid;}

class Mesh
{
public:
    Mesh(std::string filename);
    ~Mesh();

    void writeMesh(std::string filename, std::string value_label, const std::vector<double>& values) const;
    int getElement3DCount() const {return elements_3d_.size(); }

    void partitionMesh();                           //write this function
	void populateMeshVertices(double alive_probability);
	void updateVertexStates();
    void outputStatistics() const;                  //write this function
	void writeCellMatrixFile(
		std::string filename,
		bool append
	);
	void calculateVertexConnectivity();
	double calculateVertexLifeRate();

private:

    void readMesh(std::string filename);
    void getElementVertices();

    const vertex_t* getVertexFromMID(int vertex_mid) const;
    void createElementCentroidsList();



    int rank_;
    int nproc_;
	int n_shared_vertices_;
	std::map<int, int> mid2lindx_;
    std::vector<vertex_t> vertices_;
    std::vector<element_t> elements_3d_;
    std::vector<double3_t> element_3d_centroids_;
};

#endif
