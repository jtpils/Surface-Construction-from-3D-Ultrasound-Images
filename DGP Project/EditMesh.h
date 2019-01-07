/* Copyright (c) Darcy Harisson, Russell Gillette
 * April 2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
 * to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#pragma once

#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include "ImageDataLoader.h"
//#define USE_PREV

const std::size_t HOLE_INDEX = static_cast<std::size_t>( -1 );

struct half_edge{
	std::size_t next; // Index of next half-edge in the loop.
#ifdef USE_PREV
	std::size_t prev; // Index of the previous half-edge in the loop. Alternatively: use next->next if you have strictly triangle meshes.
#endif
	std::size_t twin; // Index of half-edge that is in other face that shares this edge. Open edges?
	std::size_t vert; // Index of vertex at the start of this edge.
	std::size_t face; // Index of face to the "left" of this edge.
};

class EditMesh;

/* an iterator over the one ring of vertices pointed to by this
 * half edge 
 */
class vvert_iterator{
private:
	friend class EditMesh;
	const half_edge* m_cur;
	const half_edge* m_end;
};

class vface_iterator{
//private:
public:
	friend class EditMesh;
	const half_edge* m_cur;
	const half_edge* m_end;
	const half_edge* m_next;
};

class fvert_iterator{
private:
	friend class EditMesh;
	const half_edge* m_cur;
	const half_edge* m_end;
};

class fface_iterator{
private:
	friend class EditMesh;
	const half_edge* m_cur;
	const half_edge* m_end;
};
struct vertErr
{
  public:
    double _vertError;
    size_t _vertI;
    vertErr( double vertError, int vertI ) : _vertError(vertError), _vertI(vertI) {}
};
struct userConstraint
{
  public:
    Eigen::Vector3d _vertPos;
    size_t _vertI;
    userConstraint( size_t vertI, Eigen::Vector3d vertPos ) : _vertPos(vertPos), _vertI(vertI) {}
};

struct vertices_error{
	size_t vertexIndex;
	float error;
};

struct grid_cell{
	Eigen::Vector3d pts[8];
	Eigen::Vector3d vertices[12];
	size_t vertices_global[12];
	size_t val[8];

	int cellIndex;
	int cubeindex;
	int edgeTableVal;
	int triTableVal[16];
};

struct grid{
	std::vector<size_t> row;
	std::vector<size_t> col;
	std::vector<size_t> slice;
	float isolevel;

	std::vector<grid_cell> gridCells;
	grid_cell ***gridCellMat;
};

class EditMesh{
public:
	EditMesh();

	/**
	 * Initialize the mesh from existing data.
	 * \param xyzPositions A list of doubles, storing the vertex data interleaved (ie. X1,Y1,Z1,X2,Y2,Z2,etc.)
	 * \param triangleVerts A list of vertex indices, where each run of 3 defines a triangle in the mesh. (ie. T1.v1, T1.v2, T1.v3, T2.v1, T2.v2, T2.v3, etc.)
	 */
	void init( const std::vector<double>& xyzPositions, const std::vector<std::size_t>& triangleVerts );
	void clear();

	std::size_t add_vertex( double x, double y, double z );
	std::size_t add_face( std::size_t v1, std::size_t v2, std::size_t v3 );
	std::size_t add_face( std::size_t (&f)[3] );
	
	void delete_face( std::size_t f );
	std::size_t collapse_edge( std::size_t v1, std::size_t v2 );
	
	std::size_t split_face_center( std::size_t f, std::size_t (*pOutFaceIndices)[3] = NULL );

    /**
     * returns the position of a vertex
     * \param i the vertex index (not the he index)
     * \return Vector3d the vertex position values
     */
	Eigen::Vector3d get_vertex( std::size_t i ) const;
    Eigen::Vector3d get_vnormal( std::size_t i ) const;
	//he
    Eigen::Vector3d get_fnormal( std::size_t i ) const;


	void set_vertex( std::size_t i, const Eigen::Vector3d& v );
	
    /************************************
     * Mesh modification functions
     ************************************/
public:
    /********************************
     * CS524_INPUT FUNCTIONALITY GOES HERE
     *
     * WARNING! WHEN USING EIGEN MAKE SURE
     * you check the results of the solve
     * by simply re-multiplying the solution
     * and your matrix to make sure it equals
     * your desired solution.
     *
     * There have been past Eigen bugs
     * where this did not always hold.
     ********************************/

	void Subdivision(std::size_t subdivMethod = 0);	// subdivMethod: 0 for butterfly, 1 for loop; isBoundary: 0 for closed, 1 for open mesh.
	bool isBoundaryEdge( std::size_t i ) const;
	bool isBoundaryVert( std::size_t i ) const;
	bool isBoundaryFace( std::vector< half_edge >& tempHeData, std::vector< std::size_t >& tempFaceData, std::size_t i ) const;


	//Assignment 2
	void num_vert_remove(std::size_t num = 1);
	void Errormetric_simplify(bool mode = true);
	void oneStep_simplify();
	void mesh_backward_simplify();
	size_t get_vert_number();
	size_t get_face_number();
	void get_color_vert(std::vector<size_t>& _vertVecColor);

	//the approach that didn't work
	std::vector<vertErr> UnsortedErr;
	void simplificationScheme();
	inline std::vector<size_t> EditMesh::IHeInRing(size_t vertI);
	inline std::vector<size_t> EditMesh::borderHes(size_t vertI);
	inline std::vector<half_edge> EditMesh::InnerBorderHes(size_t vertI);
	inline void EditMesh::deleteUmbrellaFaces(size_t vertI);
	bool vertRemoval(size_t vToDelete);

	// for deformation
	void get_regionOfInf_vert(std::vector<size_t>& _selected_reg_influence_verts);
	void clear_regionOfInf_vert();
	void get_control_vert(std::size_t& _control_vert);
	void clear_control_vert();
	void clear_displacementQuanta();
	void move_control_vert();

	//main algorithm - laplacian solver
	Eigen::VectorXd bX;
	Eigen::VectorXd bY;
	Eigen::VectorXd bZ;

	void popBmat(std::vector<size_t> ROI_verts, std::vector<bool> isBoundaryRegion,std::vector<size_t> ITL);
	std::vector<Eigen::Triplet<double>> popAmat(std::vector<size_t> ROI_verts, std::vector<bool> isBoundaryRegion,
										std::vector<size_t> ITL, std::vector<size_t> ITR);
	// for deformation
	//deformation algorithm

	std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > m_vertices_new;
	std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > m_vertices_init;
	
	std::vector<size_t> find_borderVerts(size_t vertI);
	
	size_t softConst;
	//Eigen::SparseMatrix<double>* A;
	Eigen::VectorXd B_x, B_y, B_z;
	//Eigen::VectorXd p_x,p_y,p_z;

	//given constraints by the user
	std::vector<userConstraint> userConstraints;
	void set_userConstraints(std::vector<size_t> _selected_reg_influence_verts);

	//columns of A are euql to the number of verts
	//rows are euql to the number of verts + user defined constraints
	void deformationScheme();
	Eigen::SparseMatrix<double> compAmat();
	void compBmat();
	size_t get_rows();
	double compWeights(std::size_t vFrom, std::size_t vTo );
	
	Eigen::MatrixXd get_Arot(std::size_t v,std::size_t vj);
	Eigen::MatrixXd get_R(std::size_t v);
	void oneStep_iteration();

	//-------------------- Simplification - Yasmin
	std::vector<Eigen::Vector3d> getNeighborsXYZ(size_t vertIndex);
	std::vector<int> getNeighborsIndices(size_t vertIndex);
	float computeErrorAtVertex(size_t vertIndex);
	void sortVerticesByError(); // sorts in descending order
	void distributeErrorToNeighbors(size_t vertIndex);
	std::vector<half_edge> getConnectedHalfEdges(size_t vertIndex); // returns the half edges pointing at the vertex, originating from the neighbors
	std::vector<size_t> getConnectedHalfEdgesIndices(size_t vertIndex);
	std::vector<size_t> getBorderVertices(size_t vertIndex);
	std::vector<half_edge> getBorderHalfEdges(std::vector<size_t> borderingVertices);
	std::vector<size_t> EditMesh::getBorderHalfEdgesIndices(std::vector<size_t> borderingVertices);
	void triangulateHole(size_t vertIndex, std::vector<size_t> borderVertices);
	void printVerticesWithErrors();
	void deleteUmbrellaFacesY(size_t vertIndex);
	void deleteFacesAndTriangulate(size_t vi);
	void computeInitialErrors_simpleMethod();
	void computeInitialErrors_advancedMethod();
	size_t getLeastErrorVertex();
	void simplify(size_t numOfVertsToRemove);
	bool isVertexRemovalValid(size_t vi);
	void moveVertextToEnd(size_t vi);
	bool isSimplest();
	size_t getNumRemainingVerts();

	std::vector<vertices_error> vertices;
	size_t numOfRemovedVert;

	/*================================================
	* mesh generator functions - Marching Cubes
	*================================================*/
	std::vector<size_t> triangles;
	std::vector<Eigen::Vector3d> verts;
	std::vector<double> vertsInput;
	double scaleFactor;
	//double scaleFactor;
	double zScaleFactor;
	//double zScaleFactor;

	grid mainGrid;
	int edgeTable[256];
	int triTable[256][16];


	void createImagesGrid(std::vector<Eigen::MatrixXd> images);
	void createGridCells(std::vector<Eigen::MatrixXd> images);
	void initializeTables();
	size_t getVertexIndex(Eigen::Vector3d vertex, size_t cubeRow, size_t cubeCol, size_t cubeSlice);
	size_t getVertexIndex(Eigen::Vector3d vertex);
	size_t getGlobalVertexIndex(Eigen::Vector3d vertex, size_t cubeCount);
	grid_cell getCell(size_t row, size_t col, size_t slice, size_t &cellIndex);
	void writeTrianglesAndVertsToFile();
	void readTrianglesAndVertsFromFileAndBuildMesh56();
	void readTrianglesAndVertsFromFileAndBuildMesh30();
	void readTrianglesAndVertsFromFileAndBuildMesh113();

	//--------------------FINAL PROJECT
	void rm_flip_edge();
	size_t newVertsNum;

	void refine_mesh();
	void smooth_mesh();
	void re_mesh();
	//difference between the moved vertex normal and the 1-ring face normals
	bool get_smoothnessError(size_t vertI);
	std::vector<Eigen::Vector3d> oldVertices_normals;

	void loadMesh();
	void load30Slices();
	void load56Slices();
	void load113Slices();
	//std::vector<size_t> triangles;
	//std::vector<Eigen::Vector3d> verts;
	//std::vector<double> vertsInput;
	//void writeTrianglesAndVertsToFile();
	//void readTrianglesAndVertsFromFileAndBuildMesh();
	/*====================================================
    * heap functions
    *====================================================*/
	inline int EditMesh::Left(int iIndex);

	inline int EditMesh::Right(int iIndex); 

	inline int EditMesh::Parent(int iIndex) ;

	inline void EditMesh::Swap(vertErr &irX, vertErr &irY) ;

	inline int EditMesh::SwapWithChild(int iIndex, vertErr* ipHeap, int iSize) ;

	inline void EditMesh::RemoveRoot(vertErr* ipHeap, int iSize) ;

	inline int EditMesh::SwapWithParent(int i, vertErr* ipHeap) ;
	inline void EditMesh::AddElement(vertErr iNewEntry, vertErr* ipHeap, int iSize) ;
	inline void EditMesh::OutputArray(vertErr* ipArray, int iBar, int iSize) ;

private:
	std::size_t m_numVert;
	std::size_t m_numColored;
	size_t m_numVertRemoved;
	std::vector<double> m_simplifyError;
	bool m_simplifyMetric;	// true for angle, false for plane
	std::vector<size_t> vertVecColor;
	std::stack<std::vector<size_t>> removed_verts;

public:	
	// for deformation
	bool IncludeRot;
	std::set<size_t> selected_reg_influence_verts;
	std::size_t control_vert;
	Eigen::Vector3d DisplacementQuanta;
	std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > m_vertices_old;

	std::vector<size_t> boundaryColor;	// for debug

private:
	bool remove_vetex(std::size_t vertIndex);
	double get_weight_sum(std::size_t vertex);
	Eigen::Vector3d get_old_delta(std::size_t vertex);
	Eigen::Vector3d get_rotated_delta(std::size_t vertex);

	
	bool remove_one_vertex();
	bool vetex_removal_algorithm(std::size_t vertIndex);
	double angle_two_vector(Eigen::Vector3d v1, Eigen::Vector3d v2 );
	double get_Gauss_error(std::size_t vertex );
	double get_angle_error(std::size_t vertex);

    /*================================================
     * Iterator functions
     *================================================*/
public:
	/**
	 * Initialize an iterator that visits the 1-ring of 'vertex'.
	 * \param it The iterator to initialize.
	 * \param vertex The index of the vertex to iterate around.
	 * \return False if the vertex has no neighbors (ie. a floating vertex)
	 */
	bool init_iterator( vvert_iterator& it, std::size_t vertex ) const;
	
    /**
     * move the iterator to a boundary, if one exists, and set its end
     * to the current location
     * \param it The iterator to move
     * \return False if there is no boundary
     */
    bool reset_boundary_iterator( vvert_iterator& it ) const;

    /**
     * \param face index
     * \return true if on boundary false else
     */
    bool isBoundaryFace( std::size_t i ) const;


	/**
	 * Advances a vertex iterator to the next vertex in the 1-ring.
	 * \param it The iterator to advance
	 * \return False if the iterator has completed the loop. It may continue to be used at this point as it will merely restart the loop.
	 */
	bool advance_iterator( vvert_iterator& it ) const;
	
	std::size_t deref_iterator( const vvert_iterator& it ) const;
	std::size_t deref_iterator_left_face( const vvert_iterator& it ) const;
	std::size_t deref_iterator_right_face( const vvert_iterator& it ) const;
	std::size_t deref_iterator_left_edge( const vvert_iterator& it ) const;
	std::size_t deref_iterator_right_edge( const vvert_iterator& it ) const;

	Eigen::Vector3d get_vertex( const vvert_iterator& it ) const;
	double get_cotan_weight( const vvert_iterator& it ) const;
	double get_weight( const vvert_iterator& it ) const;

	/**
	 * Initialize an iterator that visits the faces in the 1-ring of 'vertex'
	 */
	bool init_iterator(vface_iterator& it, std::size_t vertex ) const;
	bool advance_iterator( vface_iterator& it ) const;
	std::size_t deref_iterator( const vface_iterator& it ) const;
	Eigen::Vector3d get_normal( const vface_iterator& it ) const;
	Eigen::Vector4d get_plane( const vface_iterator& it ) const;

	bool init_iterator(fvert_iterator& it, std::size_t face ) const;
	bool advance_iterator( fvert_iterator& it ) const;
	std::size_t deref_iterator( const fvert_iterator& it ) const;

	bool init_iterator(fface_iterator& it, std::size_t face ) const;
	bool advance_iterator( fface_iterator& it ) const;
	std::size_t deref_iterator( const fface_iterator& it ) const;

    /*====================================================
     * Helper Functions
     *====================================================*/
public:
    void getIndicesForFace( size_t tri_index, size_t indicesForFace[3] );
    Eigen::Vector3d getFaceMidpoint(size_t tri_index);

    /*====================================================
     * Reflection Interface (learn stuff about the mesh)
     *====================================================*/
public:
    /* get a list of all vertices/face indices in the mesh as floats in contiguous
     * memory. Used for rendering so loss of precision unimportant */
    void get_draw_data( float *verts, int *indices ) const;
    void get_draw_normals( float *normals ) const;
    int  get_edit_count() const;
    void get_face_neighbors(int face_index, size_t neighbors[3]);
	void edit_count_updater();

    /* get information about internal state */
    std::size_t get_vert_size() const;
    std::size_t get_face_size() const;

    void test_flip();
	static void test();
	void verify() const;

	void write_to_obj_stream( std::ostream& stream ) const;

public:
	const half_edge& prev( const half_edge& cur ) const;
	const half_edge& next( const half_edge& cur ) const;
	const half_edge& twin( const half_edge& cur ) const;

private:
	half_edge* find_edge( std::size_t vFrom, std::size_t vTo );
	half_edge* find_twin( std::size_t vFrom, std::size_t vTo );
	std::size_t collapse_edge( std::size_t he );

	void delete_half_edge_impl( std::size_t he );

	template <std::size_t N>
	void delete_half_edges_impl( std::size_t (&edges)[N] );
	
	// Splits a boundary edge into 3 by adding 2 vertices on the specified edge, connecting them to the other vertex. Replaces this face with 3 new ones.
	void split_boundary_edge( std::size_t he, std::size_t (*pOutVertIndices)[2] = NULL, std::size_t (*pOutFaceIndices)[3] = NULL );

    /**
     * a simple helper function to determine if adding a particular face
     * will break mesh manifoldness
     * \param input vertices that could become a face
     * \return if this face will break manifoldness
     */
    bool is_safe_addface( std::size_t v1, std::size_t v2, std::size_t v3 );
    /**
     * Flip the passed in edge to connect the two opposing vertices
     * Flips both the passed in half_edge and its twin
     * \param cur The half edge to flip
     * \return True if successfully flipped, false if edge
     */
    bool flip_edge( half_edge &cur );

    /* =========================================
     * Mesh Collision Tests
     * =========================================*/
public:
    void updateBBox();

    // TODO: update when mesh is changed
    // NOTE: as of right now only used by skeleton animation code
    Eigen::Vector3d bboxMin;
    Eigen::Vector3d bboxMax;
    Eigen::Vector3d bSphereCenter;
    double          bSphereRadius;

    /*========================================
     * Variables for Specific Functions
     *========================================*/
private:
    unsigned int subdiv_iter; // used to track even applications of subdivision

    /*****************************************
     * Mesh Variables
     *****************************************/
private:
	friend class mesh_adjacency;
    int edit_count; // increment this value every change

	// TODO: Switch to std::deque to avoid pointer invalidation when adding new half-edges.
	std::vector< half_edge > m_heData;     // All the half-edges that make up the mesh.
	std::vector< std::size_t > m_faceData; // A mapping from face index to an arbitrary half-edge on its boundary.
	std::vector< std::size_t > m_vertData; // A mapping from vertex index to an arbitrary half-edge originating from this vertex. Can be "HOLE_INDEX" for unconnected vertices.
	
	std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > m_vertices;
	std::vector< Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d> > m_simplifyData;
	
	struct node{
		double value;
	};

	std::vector<node> m_simplifyQueue;

	// TODO: If we need arbitrary data, add a std::map< std::string, vector > that holds the per-vertex data. 
	// TODO: Same for per-half-edge and per-face data.
	friend void collapse_one_edge( EditMesh& );
	bool add_back_one_vertex();

};

void collapse_one_edge( EditMesh& m );

inline EditMesh::EditMesh() 
	: edit_count(0), subdiv_iter(0), m_numVertRemoved(0), 
	m_simplifyError(true), m_numVert(1), m_numColored(10), 
	control_vert(std::numeric_limits<size_t>::max()), IncludeRot(true)
{}

inline void EditMesh::clear(){
	subdiv_iter = 0;
	edit_count = 0;

	m_heData.clear();
	m_faceData.clear();
	m_vertData.clear();
	m_vertices.clear();
}

inline std::size_t EditMesh::add_vertex( double x, double y, double z ){
	std::size_t newIndex = m_vertices.size();
	m_vertices.emplace_back( x,y,z );
	m_vertData.push_back( HOLE_INDEX );
	return newIndex;
}

inline Eigen::Vector3d EditMesh::get_vertex( std::size_t i ) const {
	return m_vertices[i];
}

inline void EditMesh::set_vertex( std::size_t i, const Eigen::Vector3d& v ) {
	m_vertices[i] = v;
}

inline std::size_t EditMesh::add_face( std::size_t v1, std::size_t v2, std::size_t v3 ) {
	std::size_t v[] = { v1, v2, v3 };
	return this->add_face( v );
}

inline std::size_t EditMesh::collapse_edge( std::size_t v1, std::size_t v2 ){
	return this->collapse_edge( this->find_twin( v1, v2 )->twin );
}

inline bool EditMesh::init_iterator( vvert_iterator& it, std::size_t vertex ) const {
	std::size_t vertToHE = m_vertData[ vertex ];
			
	if( vertToHE == HOLE_INDEX ) // float point
		return false;

	// Store a pointer to the (arbitrary) first half-edge pointing into the specified vertex. This implies m_heData[m_cur->twin].vert == vertex & m_heData[m_cur->next].vert == vertex.
	// By iterating around the half-edges pointing into 'vertex' we will visit all the vertices in the 1-ring.
	it.m_cur = it.m_end = &m_heData[ m_heData[ vertToHE ].twin ];
	return true;
}

inline bool EditMesh::reset_boundary_iterator( vvert_iterator &it ) const {
   const half_edge *cur = it.m_cur;
    do {
        if (it.m_cur->face == HOLE_INDEX)
        {
            it.m_end = it.m_cur;
            return true;
        }
        it.m_cur = &m_heData[ m_heData[ it.m_cur->next ].twin ];
    } while (cur != it.m_cur);

    return false;
}

inline bool EditMesh::isBoundaryFace( std::size_t i ) const {
    std::size_t he_index = m_faceData[i];
    const half_edge *he = &m_heData[ he_index ];
    while( he->next != he_index ) {
        if( m_heData[he->twin].face == HOLE_INDEX )
            return true;
        he = &m_heData[ he->next ];
    }
    return false;
}

inline bool EditMesh::advance_iterator( vvert_iterator& it ) const {
	it.m_cur = &m_heData[m_heData[ it.m_cur->next].twin];
	return it.m_cur != it.m_end;
}

inline std::size_t EditMesh::deref_iterator( const vvert_iterator& it ) const {
	return it.m_cur->vert;
}

inline std::size_t EditMesh::deref_iterator_left_face( const vvert_iterator& it ) const {
	return m_heData[it.m_cur->twin].face;
}

inline std::size_t EditMesh::deref_iterator_right_face( const vvert_iterator& it ) const {
	return it.m_cur->face;
}

inline std::size_t EditMesh::deref_iterator_left_edge( const vvert_iterator& it ) const {
	return it.m_cur->twin;
}

inline std::size_t EditMesh::deref_iterator_right_edge( const vvert_iterator& it ) const {
	return m_heData[ it.m_cur->twin ].twin;
}

inline bool EditMesh::init_iterator( vface_iterator& it, std::size_t vertex ) const {
	std::size_t vertToHE = m_vertData[ vertex ];
			
	if( vertToHE == HOLE_INDEX )
		return false;

	// Store a pointer to the (arbitrary) first half-edge pointing into the specified vertex. This implies m_heData[m_cur->twin].vert == vertex & m_heData[m_cur->next].vert == vertex.
	// By iterating around the half-edges pointing into 'vertex' we will visit all the vertices in the 1-ring.
	it.m_cur = it.m_end = &m_heData[ m_heData[ vertToHE ].twin ];
	it.m_next = &m_heData[ m_heData[ it.m_cur->next ].twin ];
	return true;
}

inline bool EditMesh::advance_iterator( vface_iterator& it ) const {
	it.m_cur = it.m_next;
	it.m_next = &m_heData[ m_heData[ it.m_next->next ].twin ];
	return it.m_cur != it.m_end;
}

inline std::size_t EditMesh::deref_iterator( const vface_iterator& it ) const {
	return it.m_cur->face;
}

inline bool EditMesh::init_iterator(fvert_iterator& it, std::size_t face ) const {
	assert( face < m_faceData.size() );
	std::size_t faceToHE = m_faceData[ face ];
	assert( faceToHE < m_heData.size() );

	// Store a pointer to the (arbitrary) first half-edge pointing into the specified vertex. This implies m_heData[m_cur->twin].vert == vertex & m_heData[m_cur->next].vert == vertex.
	// By iterating around the half-edges pointing into 'vertex' we will visit all the vertices in the 1-ring.
	it.m_cur = it.m_end = &m_heData[ faceToHE ];
	return true;
}

inline bool EditMesh::advance_iterator( fvert_iterator& it ) const {
	it.m_cur = &m_heData[ it.m_cur->next ];
	return it.m_cur != it.m_end;
}

inline std::size_t EditMesh::deref_iterator( const fvert_iterator& it ) const {
	return it.m_cur->vert;
}

inline bool EditMesh::init_iterator(fface_iterator& it, std::size_t face ) const {
	assert( face < m_faceData.size() );
	std::size_t faceToHE = m_faceData[ face ];
	assert( faceToHE < m_heData.size() );

	// Store a pointer to the (arbitrary) first half-edge pointing into the specified vertex. This implies m_heData[m_cur->twin].vert == vertex & m_heData[m_cur->next].vert == vertex.
	// By iterating around the half-edges pointing into 'vertex' we will visit all the vertices in the 1-ring.
	it.m_cur = it.m_end = &m_heData[ faceToHE ];
	return true;
}

inline bool EditMesh::advance_iterator( fface_iterator& it ) const {
	it.m_cur = &m_heData[ it.m_cur->next ];
	return it.m_cur != it.m_end;
}

inline std::size_t EditMesh::deref_iterator( const fface_iterator& it ) const {
	return m_heData[ it.m_cur->twin ].face;
}

inline const half_edge& EditMesh::prev( const half_edge& cur ) const {
#ifdef USE_PREV
	return m_heData[ cur.prev ];
#else
	return m_heData[ m_heData[ cur.next ].next ];
#endif
}

inline const half_edge& EditMesh::next( const half_edge& cur ) const {
	return m_heData[ cur.next ];
}

inline const half_edge& EditMesh::twin( const half_edge& cur ) const {
	return m_heData[ cur.twin ];
}

inline std::size_t EditMesh::get_vert_size() const {
    return m_vertData.size();
}

inline std::size_t EditMesh::get_face_size() const {
    return m_faceData.size();
}

inline int EditMesh::get_edit_count() const{
    return edit_count;
}

inline void EditMesh::get_face_neighbors(int face_index, size_t neighbors[3]) {
    fface_iterator fit;
    init_iterator(fit, face_index);
    int i = 0;
    do {
        neighbors[i++] = deref_iterator(fit);
    } while(advance_iterator(fit));

    for (; i < 3; i++)
        neighbors[i] = HOLE_INDEX;
}

inline void EditMesh::edit_count_updater(){
	++edit_count;
}

inline bool EditMesh::isBoundaryFace( std::vector< half_edge >& tempHeData, std::vector< std::size_t >& tempFaceData, std::size_t i ) const
{
	std::size_t he_index = tempFaceData[i];
	const half_edge *he = &tempHeData[ he_index ];
	while( he->next != he_index ) {
		if( tempHeData[he->twin].face == HOLE_INDEX )
			return true;
		he = &tempHeData[ he->next ];
	}
	return false;
}

inline bool EditMesh::isBoundaryEdge( std::size_t i ) const
{
	const half_edge *he = &m_heData[i];
	if (he->face == HOLE_INDEX || m_heData[he->twin].face == HOLE_INDEX)
		return true;
	return false;
}

inline bool EditMesh::isBoundaryVert( std::size_t i ) const
{
	vvert_iterator vvIter;
	
	if (!this->init_iterator( vvIter, i))
		std::cerr << "error because of float point" << std::endl;
	
	do 
	{
		if (this->isBoundaryEdge(vvIter.m_cur->twin))
			return true;
		
	} while (this->advance_iterator(vvIter));

	return false;
}

inline size_t EditMesh::get_vert_number()
{
	return m_vertData.size();
}

inline size_t EditMesh::get_face_number()
{
	return m_faceData.size();
}
inline int EditMesh::Left(int iIndex) {
	return ((iIndex << 1) + 1);
}

inline int EditMesh::Right(int iIndex) {
	return ((iIndex << 1) + 2);
}

inline int EditMesh::Parent(int iIndex) {
	return ((iIndex - 1) >> 1);
}

inline void EditMesh::Swap(vertErr& irX, vertErr& irY) {
	vertErr iTemp = irX;
	irX = irY;
	irY = iTemp;
}

inline int EditMesh:: SwapWithChild(int iIndex, vertErr* ipHeap, int iSize) {
	int iLeft		= Left(iIndex);
	int iRight		= Right(iIndex);
	int iLargest	= iIndex;
	if (iRight < iSize) {
		if (ipHeap[iLeft]._vertError < ipHeap[iRight]._vertError) {
			iLargest = iRight;
		} else {
			iLargest = iLeft;
		}
		if (ipHeap[iIndex]._vertError > ipHeap[iLargest]._vertError) {
			iLargest = iIndex;
		}
	} else if (iLeft < iSize) {
		if (ipHeap[iIndex]._vertError < ipHeap[iLeft]._vertError) {
			iLargest = iLeft;
		}
	}
	if (ipHeap[iIndex]._vertError < ipHeap[iLargest]._vertError) {
		Swap(ipHeap[iIndex], ipHeap[iLargest]);
	}
	return iLargest;
}

inline void EditMesh::RemoveRoot(vertErr* ipHeap, int iSize) {
	// Swap the last element with the root
	Swap(ipHeap[0], ipHeap[iSize - 1]);
	--iSize;
	int iLasti = 0;
	int i = SwapWithChild(0, ipHeap, iSize);
	while (i != iLasti) {
		iLasti = i;
		i = SwapWithChild(i, ipHeap, iSize);
	}
}

inline int EditMesh::SwapWithParent(int i, vertErr* ipHeap) {
	if (i < 1) {
		return i;
	}
	int iParent = Parent(i);
	if (ipHeap[i]._vertError > ipHeap[iParent]._vertError) {
		Swap(ipHeap[i], ipHeap[iParent]);
		return iParent;
	} else {
		return i;
	}
}
inline void EditMesh::AddElement(vertErr iNewEntry, vertErr* ipHeap, int iSize) {
	ipHeap[iSize] = iNewEntry;
	int iLasti = iSize;
	int i = SwapWithParent(iLasti, ipHeap);
	while (iLasti != i) {
		iLasti = i;
		i = SwapWithParent(i, ipHeap);
	}
}


inline void EditMesh::OutputArray(vertErr* ipArray, int iBar, int iSize) {
	using namespace std;
	for (int i = 0; i < iSize; ++i) {
		if (i == iBar) {
			//cout << "|  ";
		}
		//cout << ipArray[i]._vertError << "  ";
	}
	//cout << endl;
}

//Store half edges originating form  the specified vertex.
//CCW
inline std::vector<size_t> EditMesh::IHeInRing(size_t vertI){
	vvert_iterator it;
	std::vector< size_t > heVec;
		if( this->init_iterator( it, vertI ) ){

			std::size_t HeInRingInd = it.m_cur->twin;
			half_edge HeInRing		= m_heData[it.m_cur->twin];
			std::size_t firstHeInRingInd = HeInRingInd;
			
			//std::cout<<"firstHeInRingInd is: "<<firstHeInRingInd<<std::endl;
			do{
				//printf("HeInRingInd is %d \n", HeInRingInd);
				heVec.push_back(HeInRingInd);
				//std::cout<<HeInRingInd;
				//std::cout<<std::endl;

				HeInRingInd =   m_heData[m_heData[HeInRing.next].next].twin;
				
				HeInRing	=  m_heData[HeInRingInd];
				if (heVec.size()>7){
					std::cout<<"broke here"<<std::endl;
					break;}
			}while( firstHeInRingInd  != HeInRingInd);
			
			//printf("There are %d edges in 1ring of the vertex\n", heVec.size()-1);
			return heVec;
		}
}
inline std::vector<size_t> EditMesh::borderHes(size_t vertI){
	std::vector<size_t>heBorder;
	size_t k = 1;
	size_t heBorderNext = m_heData[m_heData[m_vertData[vertI]].next].twin;
	do{
		heBorder.push_back(heBorderNext);
		heBorderNext = m_heData[m_heData[m_heData[m_heData[m_heData[m_heData[ m_heData[ heBorder[k-1] ]
		.twin ].next].next].twin].next].next].twin;		
		k++;
	}while( heBorderNext != heBorder[0]);
	/*std::cout<<"Number of border half-edges originating from this vert is"<<" "<<heBorder.size();
	std::cout<<std::endl;*/
	return heBorder;
}

inline std::vector<half_edge> EditMesh::InnerBorderHes(size_t vertI){
	std::vector<half_edge>heBorder;
	size_t k = 0;
	size_t heBorderNext = m_heData[m_vertData[vertI]].next;
	size_t iHeBorder0	= heBorderNext;
	//std::cout<<"inner border hes twins are: "<<std::endl;
	do{
		heBorder.push_back( m_heData[heBorderNext]);
		heBorderNext = m_heData[m_heData[ m_heData[ heBorderNext ]
		.next].twin].next;		
		//std::cout<<heBorder[k].twin<<std::endl;
		k++;
	}while( heBorderNext != iHeBorder0);
	//std::cout<<"Number of border half-edges originating from this vert is"<<" "<<heBorder.size();
	//std::cout<<std::endl;
	return heBorder;
}
inline void EditMesh::deleteUmbrellaFaces(size_t vertI){
	size_t he;
	while(m_vertData[vertI]!= HOLE_INDEX){
		he = m_vertData[vertI];
		delete_face(m_heData[he].face);
		
	}
}

inline std::vector<size_t> EditMesh::find_borderVerts(size_t vertI){
	// initialize the iterator visit the one ring of "vertex"
	vvert_iterator it;
	if (!this->init_iterator(it, vertI))
		std::cerr << "error because of float point" << std::endl;
	std::vector<size_t> borderVerts;
	// sum the 1-ring angles	
	do {
		borderVerts.push_back(it.m_cur->vert);

	} while(advance_iterator(it));
	return borderVerts;
}