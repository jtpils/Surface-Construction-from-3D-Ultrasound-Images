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
#define _SCL_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "EditMesh.h"
#include "Utils.h"


#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <cmath>

#include <map>
#include <set>
#include <stack>
#include <iostream>
#include <memory>
#include <fstream>
#include <algorithm>


using namespace std;


//const double PI(3.141592653589793238);

//#if defined(NDEBUG) && defined(ALWAYS_ASSERT)
//#undef NDEBUG
//#endif
//#include <cassert>

EditMesh *loadEditMeshFromFile(std::string file_name);

namespace detail{
	inline void init( half_edge& he, std::size_t next, std::size_t twin, std::size_t vert, std::size_t face ){
		he.next = next;
		he.twin = twin;
		he.vert = vert;
		he.face = face;
	}

	void delete_face( std::vector<std::size_t>& faceData, std::vector<half_edge>& heData, std::size_t f ){
		assert( f < faceData.size() );

		// In order to delete the face properly, we need to move a face from the end of the list to overwrite 'f'. Then we need to update the 
		// indices stored in the moved face's half-edges.
		faceData[f] = faceData.back();
		faceData.pop_back();

		if( f != faceData.size() ){
			//std::clog << "Reindexed face " << faceData.size() << " to " << f << std::endl;

			std::size_t he = faceData[f];
			do {
				assert( heData[he].face == faceData.size() );

				heData[he].face = f;
				he = heData[he].next;
			} while( he != faceData[f] );
		}
	}

	template <int N>
	void delete_faces( std::vector<std::size_t>& faceData, std::vector<half_edge>& heData, std::size_t (&fToDelete)[N] )
	{
		// Sort the faces by decreasing index so that we can safely delete them all without causing any of them to be accidentally re-indexed (which
		// cause 'fToDelete' to contain invalid indices). This also chooses the optimal deletion order to minimize re-indexing.
		std::sort( fToDelete, fToDelete + N, std::greater<std::size_t>() );
		for( std::size_t i = 0; i < N; ++i )
			detail::delete_face( faceData, heData, fToDelete[i] );
	}
}

void init_adjacency( std::size_t numVertices, const std::vector<std::size_t>& faces, std::vector< half_edge >& m_heData, std::vector< std::size_t >& m_faceData, std::vector< std::size_t >& m_vertData ){
	typedef std::map< std::pair<std::size_t, std::size_t>, std::size_t > edge_map_type;
	
	assert( faces.size() % 3 == 0 && "Invalid data specified for faces. Must have 3 vertex indices per face." );

	edge_map_type edgeMap; // Use a temporary map to find edge pairs.

	m_heData.reserve( faces.size() ); // Assume there are 3 edges per face.
	m_faceData.resize( faces.size() / 3 );
	m_vertData.resize( numVertices, HOLE_INDEX ); // Init with HOLE_INDEX since a vert might be floating w/ no fac   es.

	for( std::size_t i = 0, iEnd = faces.size(); i < iEnd; i+=3 ){
		std::size_t f[] = { faces[i], faces[i+1], faces[i+2] };
		std::size_t fIndex = i / 3;

		// The index of the first (of three) half-edges associated with the current face.
		std::size_t heIndex = m_heData.size();

		half_edge he[3];
		detail::init( he[0], heIndex+1, HOLE_INDEX, f[0], fIndex );
		detail::init( he[1], heIndex+2, HOLE_INDEX, f[1], fIndex );
		detail::init( he[2], heIndex, HOLE_INDEX, f[2], fIndex );
#ifdef USE_PREV
		he[0].prev = heIndex+2;
		he[1].prev = heIndex;
		he[2].prev = heIndex+1;
#endif
			
		// These will be set each time a vertex is referenced, but that's fine. The last assignment will stick.
		m_faceData[ fIndex ] = heIndex;
		m_vertData[ f[0] ] = heIndex;
		m_vertData[ f[1] ] = heIndex+1;
		m_vertData[ f[2] ] = heIndex+2;

		edge_map_type::iterator it;

		it = edgeMap.lower_bound( std::make_pair( f[0], f[1] ) );
		if( it != edgeMap.end() && it->first.first == f[0] && it->first.second == f[1] ){
			m_heData[it->second].twin = heIndex;
			he[0].twin = it->second;
			edgeMap.erase( it );
		}else{
			he[0].twin = HOLE_INDEX;
			edgeMap.insert( it, std::make_pair( std::make_pair( f[1], f[0] ), heIndex ) ); // NOTE: Reversed order since we are matching opposite half_edge.
		}

		it = edgeMap.lower_bound( std::make_pair( f[1], f[2] ) );
		if( it != edgeMap.end() && it->first.first == f[1] && it->first.second == f[2] ){
			m_heData[it->second].twin = heIndex+1;
			he[1].twin = it->second;
			edgeMap.erase( it );
		}else{
			he[1].twin = HOLE_INDEX;
			edgeMap.insert( it, std::make_pair( std::make_pair( f[2], f[1] ), heIndex+1 ) ); // NOTE: Reversed order since we are matching opposite half_edge.
		}

		it = edgeMap.lower_bound( std::make_pair( f[2], f[0] ) );
		if( it != edgeMap.end() && it->first.first == f[2] && it->first.second == f[0] ){
			m_heData[it->second].twin = heIndex+2;
			he[2].twin = it->second;
			edgeMap.erase( it );
		}else{
			he[2].twin = HOLE_INDEX;
			edgeMap.insert( it, std::make_pair( std::make_pair( f[0], f[2] ), heIndex+2 ) ); // NOTE: Reversed order since we are matching opposite half_edge.
		}

		m_heData.push_back( he[0] );
		m_heData.push_back( he[1] );
		m_heData.push_back( he[2] );
	}

	// Keep track of the last edge we processed so we can hook up half_edge::prev as we go.
	std::size_t prev = HOLE_INDEX;

	// Add half-edges for any holes. Any edges still in the map are holes.
	edge_map_type::iterator it = edgeMap.begin();
	while( it != edgeMap.end() ){
		half_edge he;
		detail::init( he, HOLE_INDEX, it->second, it->first.first, HOLE_INDEX );
#ifdef USE_PREV
		he.prev = prev;
		prev = m_heData.size(); // Size is the index of the half_edge we are about to push into the list.
#endif

		m_heData[he.twin].twin = m_heData.size();
		m_heData.push_back( he );

		std::size_t curVert = it->first.first;
		std::size_t nextVert = it->first.second; // We are about to erase this information, so store it to use later.

		edgeMap.erase( it ); // We are done with this edge now.

		half_edge* twinPrev = &m_heData[m_heData[m_heData[he.twin].next].next];
		while( twinPrev->twin != HOLE_INDEX && m_heData[twinPrev->twin].face != HOLE_INDEX ){
			assert( m_heData[twinPrev->next].vert == nextVert );
			assert( m_heData[twinPrev->twin].vert == nextVert );
			twinPrev = &m_heData[m_heData[m_heData[twinPrev->twin].next].next];
		}

		if( twinPrev->twin == HOLE_INDEX ){
			// We haven't processed the next edge in the loop yet. Let's do so now so we can assume the index of the next half-edge.
			m_heData.back().next = m_heData.size();
			it = edgeMap.find( std::make_pair( nextVert, twinPrev->vert ) );
				
			assert( it != edgeMap.end() );
		}else{
			assert( m_heData[twinPrev->twin].vert == nextVert );
			assert( m_heData[twinPrev->twin].face == HOLE_INDEX );

			// We already processed this edge and have a valid index for the next half_edge.
			m_heData.back().next = twinPrev->twin;
#ifdef USE_PREV
			m_heData[ twinPrev->twin ].prev = prev; // Complete the loop
			prev = HOLE_INDEX;
#endif
			it = edgeMap.begin(); // Arbitrarily pick the next edge in the list.
		}
	}

	assert( edgeMap.empty() );
}

void EditMesh::init( const std::vector<double>& xyzPositions, const std::vector<std::size_t>& triangleVerts ){
	assert( xyzPositions.size() % 3 == 0 && "Invalid vertex positions for EditMesh::init(). Must have 3 values per-vertex." );
	assert( triangleVerts.size() % 3 == 0 && "Invalid face data for EditMesh::init(). Must have 3 vertex indices per face." );

	//// clear
	this->clear();

	//m_vertices.resize( Eigen::NoChange, xyzPositions.size() / 3 );
	m_vertices.resize( xyzPositions.size() / 3 );

	// The Eigen matrix has the same format as the incoming vector so we can straight copy it.
	// HACK: This is pretty sketchy and relies on Eigen::Vector3d having the same layout as a double[3] and nothing extra or fancy alignment.
	std::copy( xyzPositions.begin(), xyzPositions.end(), m_vertices.front().data() );

	init_adjacency( xyzPositions.size() / 3, triangleVerts, m_heData, m_faceData, m_vertData );

}

half_edge* EditMesh::find_twin( std::size_t vFrom, std::size_t vTo ){
	vvert_iterator it;
	if( !this->init_iterator( it, vFrom ) )
		return NULL;
	
	do{
		if( this->deref_iterator( it ) == vTo )
			return const_cast<half_edge*>( it.m_cur ); // Gross. This is just laziness.
	}while( this->advance_iterator( it ) );

	return NULL;
}

half_edge* EditMesh::find_edge( std::size_t vFrom, std::size_t vTo ){
	if( half_edge* he = this->find_twin( vFrom, vTo ) )
		return &m_heData[ he->twin ];
	return NULL;
}

bool EditMesh::flip_edge( half_edge &he ){
    half_edge &twin = m_heData[ he.twin ];

    if (he.face == HOLE_INDEX ||
        twin.face == HOLE_INDEX)
        return false;

    std::size_t he_tri[3];
    std::size_t twin_tri[3];

    // prep: gather half edge indices in
    // the order they should be after flip
    he_tri[0]   = he.next;
    twin_tri[0] = twin.next;
    he_tri[1]   = twin.twin;
    twin_tri[1] = he.twin;
    he_tri[2]   = m_heData[ twin_tri[0] ].next;
    twin_tri[2] = m_heData[ he_tri[0] ].next;

	if( m_heData[ he_tri[2] ].vert == m_heData[ twin_tri[2] ].vert )
		return false;

    // step 1: ensure he's verts don't point to
    // either half_edge (does not break mesh)
    m_vertData[ he.vert ] = twin_tri[0];
    m_vertData[ twin.vert ] = he_tri[0];

    // step 2: set the he's vert to new originating vert
    he.vert = m_heData[ twin_tri[2] ].vert;
    twin.vert = m_heData[ he_tri[2] ].vert;
    
    // step 3: ensure the faces point to one
    // of the half edges connected to them
    m_faceData[ he.face ] = he_tri[0];
    m_faceData[ twin.face ] = twin_tri[0];

    // step 4: fix two edges that will point
    // to the wrong face
    m_heData[he_tri[2]].face = he.face;
    m_heData[twin_tri[2]].face = twin.face;

    // step 5: ensure half edges point to
    // each other
    for( int i=0; i<3; ++i ) {
        m_heData[ he_tri[i] ].next = he_tri[(i+1)%3];
        m_heData[ twin_tri[i] ].next = twin_tri[(i+1)%3];
    }

    return true;
}

// IMPORTANT: Given a collection of half-edges to delete (ex. When removing a face we need to kill 2, 4, or 6 half-edges) they must be deleting in decreasing index order!
void EditMesh::delete_half_edge_impl( std::size_t he ){
	assert( (m_heData[he].vert >= m_vertData.size() || m_vertData[m_heData[he].vert] != he) && "Deleting this half_edge leaves a dangling link from a vertex. Must handle this first" );
		
	// Move a half_edge from the end overtop of the half_edge we are deleting, then update the indices of linked half_edges.
	m_heData[he] = m_heData.back();
	m_heData.pop_back();

	// We may have just deleted the item at the end of the list, so we have nothing to update since the indices didn't change.
	if( he != m_heData.size() )
	{
		const half_edge& heMoved = m_heData[he];

		// If the moved half_edge was the arbitrary half_edge linked to the vertex, update it.
		if( m_vertData[heMoved.vert] == m_heData.size() )
			m_vertData[heMoved.vert] = he;

		// If the moved half_edge was the arbitrary half_edge linked to the face, update it.
		if( heMoved.face != HOLE_INDEX && m_faceData[heMoved.face] == m_heData.size() )
			m_faceData[heMoved.face] = he;

		assert( heMoved.twin < m_heData.size() );
		assert( m_heData[heMoved.twin].twin == m_heData.size() );
		m_heData[heMoved.twin].twin = he;

		// NOTE: If we are deleting a bundle of half_edges, then by definition we must call delete_half_edge() in decreasing order of indices. That prevents
		//       me from having to worry about moving a partially destroyed half_edge into the 'he' position.

#ifdef USE_PREV
		assert( m_heData[heMoved.prev].next == m_heData.size() );
		m_heData[heMoved.prev].next = he;

		assert( m_heData[heMoved.next].prev == m_heData.size() );
		m_heData[heMoved.next].prev = he;
#else
		// Have to loop around the face until we find the half_edge using 'heMoved' as its 'next' entry, then update it.
		std::size_t hePrev = heMoved.next;
		while( m_heData[hePrev].next != m_heData.size() )
			hePrev = m_heData[hePrev].next;

		assert( m_heData[hePrev].next == m_heData.size() );
		m_heData[hePrev].next = he;
#endif
	}

	// Update the links in the simplification queue too.
	if( !m_simplifyQueue.empty() ){
		m_simplifyQueue[he] = m_simplifyQueue.back();
		m_simplifyQueue.pop_back();
	}
}

template <std::size_t N>
void EditMesh::delete_half_edges_impl( std::size_t (&heToDelete)[N] ){
	std::sort( heToDelete, heToDelete + N, std::greater<std::size_t>() );
	for( std::size_t i = 0; i < N; ++i )
		this->delete_half_edge_impl( heToDelete[i] );
}

bool g_debug = false;

std::size_t EditMesh::collapse_edge( std::size_t he ){
	assert( he < m_heData.size() );
	assert( m_heData[he].face != HOLE_INDEX && m_heData[m_heData[he].twin].face != HOLE_INDEX && "Cannot collapse a boundary edge" );

	const half_edge& heBase = m_heData[he];
	const half_edge& heTwin = m_heData[heBase.twin];

	// We are going to delete the faces on either side of the chosen edge, so we need to delete 3 half_edges and patch up the twin links on the 4
	// bordering edges.
	std::size_t heBorder[4];
	heBorder[0] = m_heData[ heBase.next ].twin;
	heBorder[1] = m_heData[ m_heData[ heBase.next ].next ].twin;
	heBorder[2] = m_heData[ m_heData[ heTwin.next ].next ].twin;
	heBorder[3] = m_heData[ heTwin.next ].twin;

	// TODO: Relax this assertion. We should be able to collapse a spike jutting into a hole.
	assert( ( m_heData[ heBorder[0] ].face != HOLE_INDEX || m_heData[ heBorder[1] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );
	assert( ( m_heData[ heBorder[2] ].face != HOLE_INDEX || m_heData[ heBorder[3] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );

	// Check if we can actually collapse. This checks for a degree 3 vertex not on the edge we are collapsing.
	if( m_heData[ m_heData[ m_heData[ heBorder[1] ].next ].twin ].next == heBorder[0] )
		return HOLE_INDEX;
	if( m_heData[ m_heData[ m_heData[ heBorder[2] ].next ].twin ].next == heBorder[3] )
		return HOLE_INDEX;

	// Capture the indices of things (2 faces & 6 half-edges) we want to delete.
	std::size_t fToDelete[] = { heBase.face, heTwin.face };
	std::size_t heToDelete[] = { he, heBase.next, m_heData[ heBase.next ].next, heBase.twin, heTwin.next, m_heData[ heTwin.next ].next };
	
//#ifndef NDEBUG
//	// We can't be deleting border edges!
//	for( auto i : heToDelete ){
//		if( std::find( heBorder, heBorder + 4, i ) != heBorder + 4 )
//			return HOLE_INDEX;	
//		//assert( std::find( heBorder, heBorder + 4, i ) == heBorder + 4 );
//	}
//
//	if( g_debug ){
//		std::vector< std::set<std::size_t> > verts( 3 );
//
//		verts[0].insert( heBase.vert );
//		verts[0].insert( heTwin.vert );
//
//		for( int i = 1; i < verts.size(); ++i ){
//			for( auto v : verts[i-1] ){
//				vvert_iterator it;
//				this->init_iterator( it, v );
//				do{
//					verts[i].insert( this->deref_iterator( it ) );
//				}while( this->advance_iterator( it ) );
//			}
//		}
//
//		std::vector<std::size_t> orderedVerts( verts.back().begin(), verts.back().end() );
//		std::set<std::size_t> faces;
//
//		std::vector< double > vpos;
//		std::vector< std::size_t > finds;
//
//		for( auto v : orderedVerts ){
//			vpos.push_back( m_vertices[v].x() ); vpos.push_back( m_vertices[v].y() ); vpos.push_back( m_vertices[v].z() );
//			//std::clog << "m.add_vert( " << m_vertices[v].x() << ", " << m_vertices[v].y() << ", " << m_vertices[v].z() << " );" << std::endl;
//		}
//
//		// Visit the 1-ring
//		for( auto v : verts[1] ){
//			vface_iterator it;
//			this->init_iterator( it, v );
//			do{
//				if( this->deref_iterator( it ) != HOLE_INDEX && faces.find( this->deref_iterator( it ) ) == faces.end() ){
//					faces.insert( this->deref_iterator( it ) );
//
//					fvert_iterator itFace;
//					this->init_iterator( itFace, this->deref_iterator( it ) );
//
//					std::size_t f[3];
//					std::size_t i = 0;
//					do{
//						f[i++] = std::find( orderedVerts.begin(), orderedVerts.end(), this->deref_iterator( itFace ) ) - orderedVerts.begin();
//					}while( this->advance_iterator( itFace ) );
//
//					finds.push_back( f[0] ); finds.push_back( f[1] ); finds.push_back( f[2] );
//					//std::clog << "m.add_face( " << f[0] << ", " << f[1] << ", " << f[2] << " );" << std::endl;
//				}	
//			}while( this->advance_iterator( it ) );
//		}
//
//		std::size_t base = std::find( orderedVerts.begin(), orderedVerts.end(), heBase.vert ) - orderedVerts.begin();
//		std::size_t twin = std::find( orderedVerts.begin(), orderedVerts.end(), heTwin.vert ) - orderedVerts.begin();
//		std::clog << "m.collapse_edge( " << base << ", " << twin << " );" << std::endl;
//
//		EditMesh m;
//		m.init( vpos, finds );
//		std::ofstream fout( "debug.obj" );
//		m.write_to_obj_stream( fout );
//		fout.close();
//	}
//#endif

	// We may also need to fix the vertex->half_edge link for the verts using these faces. There are technically 4, but we only update the 3 that are not going to be deleted.
	std::size_t verts[] = { this->prev( heBase ).vert, heBase.vert, this->prev( heTwin ).vert };

	// Move the base vertex (arbitrarily) to the middle of the edge. Could leave it where it is, or do something fancier too.
	//m_vertices[heBase.vert] = 0.5 * ( m_vertices[heBase.vert] + m_vertices[heTwin.vert] ); 
	m_vertices[heBase.vert] = m_vertices[heTwin.vert];

	// Adjust all the twin's 1-ring to link to the vertex we are not going to delete.
	std::size_t heIt = this->twin(this->next(heBase)).next;
	std::size_t heEnd = heBase.twin;
	for( ; heIt != heEnd; heIt = this->twin( m_heData[heIt] ).next ){
		assert( m_heData[heIt].vert == heTwin.vert );
		
		// Associate to the other vertex now, so we can delete this one.
		m_heData[heIt].vert = heBase.vert;
	}

	// Fix the vert associations if required, picking a non-hole face.
	if( m_vertData[ verts[0] ] == m_heData[ heBorder[1] ].twin )
		m_vertData[ verts[0] ] = (m_heData[ heBorder[0] ].face != HOLE_INDEX) ? heBorder[0] : m_heData[ heBorder[1] ].next;
	if( m_vertData[ verts[1] ] == he || m_vertData[ verts[1] ] == heTwin.next )
		m_vertData[ verts[1] ] = (m_heData[ heBorder[1] ].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
	if( m_vertData[ verts[2] ] == m_heData[ heBorder[2] ].twin )
		m_vertData[ verts[2] ] = (m_heData[ heBorder[3] ].face != HOLE_INDEX) ? heBorder[3] : m_heData[ heBorder[2] ].next;

	// "Delete" the other vertex
	m_vertData[heTwin.vert] = HOLE_INDEX;

	// Collapse the two triangles bordering our chosen half-edge by connecting the opposite edges together.
	m_heData[ heBorder[0] ].twin = heBorder[1];
	m_heData[ heBorder[1] ].twin = heBorder[0];
	m_heData[ heBorder[2] ].twin = heBorder[3];
	m_heData[ heBorder[3] ].twin = heBorder[2];

	// Have to delete the faces in the proper order.
	if( fToDelete[0] < fToDelete[1] )
		std::swap( fToDelete[0], fToDelete[1] );

	this->delete_half_edges_impl( heToDelete );
	detail::delete_face( m_faceData, m_heData, fToDelete[0] );
	detail::delete_face( m_faceData, m_heData, fToDelete[1] );

	return verts[1];
}

std::size_t EditMesh::add_face( std::size_t (&f)[3] ){
	std::size_t faceIndex = m_faceData.size();
	std::size_t heIndex = m_heData.size();

	// Find the half-edges on the hole face we are filling. We must either:
	//  1. Find no half-edges, if all vertices are unconnected from the mesh.
	//  3. Find one half-edge, if one of the vertices is not connected to the existing mesh.
	//  4. Find two half-edges, if we are adding a triangle inside of a polygonal hole.
	//  2. Find three half-edges, if we are filling an existing triangular hole.
	half_edge* he[] = { 
		this->find_edge( f[0], f[1] ), 
		this->find_edge( f[1], f[2] ), 
		this->find_edge( f[2], f[0] ) };
	
	// Find the first half-edge we need to modify. This is an edge 
	std::size_t base = HOLE_INDEX;
	for( std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i ){
		if( he[i] ){
			assert( he[i]->face == HOLE_INDEX && "Non-manifold mesh detected. Cannot connect to an edge which already has two incident faces (ie. One side must be a hole)" );
			if( !he[(i+2)%3] )
				base = i;
		}
	}

	if( base == HOLE_INDEX ){
		// This triangle is not connected to any others, or we completely filled a triangular hole.
		if( he[0] /*|| he[1] || he[2]*/ ){
			assert( he[0] && he[1] && he[2] );
			assert( he[0]->face == HOLE_INDEX && he[1]->face == HOLE_INDEX && he[2]->face == HOLE_INDEX );
			assert( &m_heData[ he[0]->next ] == he[1] && &m_heData[ he[1]->next ] == he[2] && &m_heData[ he[2]->next ] == he[0] );
			
			// Update the face index of the triangular hole to convert it to a face.
			he[0]->face = he[1]->face = he[2]->face = faceIndex;
			m_faceData.push_back( he[2]->next );
		}else{
			assert( !he[0] && !he[1] && !he[2] );
			assert( m_vertData[f[0]] == HOLE_INDEX && m_vertData[f[1]] == HOLE_INDEX && m_vertData[f[2]] == HOLE_INDEX && "Non-manifold mesh detected. Cannot have two hole faces incident on a vertex." );

			// Make 3 new half-edges for the triangle, and 3 new half-edges for the hole outside of the triangle.
			half_edge newHe[6];
			detail::init( newHe[0], heIndex+1, heIndex+5, f[0], faceIndex );
			detail::init( newHe[1], heIndex+2, heIndex+4, f[1], faceIndex );
			detail::init( newHe[2], heIndex  , heIndex+3, f[2], faceIndex );
			detail::init( newHe[3], heIndex+4, heIndex+2, f[0], HOLE_INDEX );
			detail::init( newHe[4], heIndex+5, heIndex+1, f[2], HOLE_INDEX );
			detail::init( newHe[5], heIndex+3, heIndex  , f[1], HOLE_INDEX );
#ifdef USE_PREV
			newHe[0].prev = heIndex+2;
			newHe[1].prev = heIndex;
			newHe[2].prev = heIndex+1;

			newHe[3].prev = heIndex+5;
			newHe[4].prev = heIndex+3;
			newHe[5].prev = heIndex+4;
#endif

			m_vertData[ f[0] ] = heIndex;
			m_vertData[ f[1] ] = heIndex+1;
			m_vertData[ f[2] ] = heIndex+2;

			m_faceData.push_back( heIndex );
			m_heData.push_back( newHe[0] );
			m_heData.push_back( newHe[1] );
			m_heData.push_back( newHe[2] );
			m_heData.push_back( newHe[3] );
			m_heData.push_back( newHe[4] );
			m_heData.push_back( newHe[5] );
		}
	}else{
		std::size_t next = (base+1)%3, prev = (base+2)%3;
		std::size_t baseIndex = static_cast<std::size_t>( he[base] - &m_heData.front() );

		assert( !he[prev] );

		if( he[next] ){
			// We have two edges to steal from the hole, and we need to add two new half-edges
			half_edge newHe[2];
			detail::init( newHe[0], baseIndex, heIndex+1, f[prev], faceIndex );
			detail::init( newHe[1], he[next]->next, heIndex, f[base], HOLE_INDEX );

#ifdef USE_PREV
			newHe[0].prev = he[base]->next;
			newHe[1].prev = he[base]->prev;

			m_heData[ he[base]->prev ].next = heIndex + 1;
			m_heData[ he[next]->next ].prev = heIndex + 1;

			he[next]->next = heIndex;
			he[base]->prev = heIndex;
#else
			// Have to find the previous half_edge in the polygonal hole so we can point it to the new half-edge in the hole.
			half_edge* hePrev = &m_heData[ he[next]->next ];
			while( &m_heData[hePrev->next] != he[base] ){
				hePrev = &m_heData[hePrev->next];
				assert( hePrev != he[next] ); // To catch weirdness.
			}
			assert( &m_heData[hePrev->next] == he[base] );
			
			hePrev->next = heIndex + 1;
			he[next]->next = heIndex;
#endif

			// Update the face indices of the half-edges to indicate they are in a triangle now.
			he[base]->face = he[next]->face = faceIndex;

			m_faceData.push_back( heIndex );
			m_heData.push_back( newHe[0] );
			m_heData.push_back( newHe[1] );
		}else{
			assert( m_vertData[ f[prev] ] == HOLE_INDEX && "Non-manifold mesh detected. Cannot have two hole faces incident on a vertex." );

			// We have one edge to steal from the hole, and we need to add four new half-edges.
			half_edge newHe[4];
			detail::init( newHe[0], baseIndex, heIndex+2, f[prev], faceIndex );
			detail::init( newHe[1], heIndex  , heIndex+3, f[next], faceIndex );
			detail::init( newHe[2], heIndex+3, heIndex  , f[base], HOLE_INDEX );
			detail::init( newHe[3], he[base]->next, heIndex+1, f[prev], HOLE_INDEX );

#ifdef USE_PREV
			newHe[0].prev = heIndex+1;
			newHe[1].prev = baseIndex;
			newHe[2].prev = he[base]->prev;
			newHe[3].prev = heIndex+2;

			m_heData[ he[base]->prev ].next = heIndex+2;
			m_heData[ he[base]->next ].prev = heIndex+3;

			he[base]->prev = heIndex;
			he[base]->next = heIndex+1;
#else
			// Have to find the previous half_edge in the polyognal hole so we can point it to the new half-edge in the hole.
			half_edge* hePrev = &m_heData[ he[base]->next ];
			while( &m_heData[hePrev->next] != he[base] ){
				hePrev = &m_heData[hePrev->next];
				assert( hePrev != he[next] ); // To catch weirdness.
			}
			assert( &m_heData[hePrev->next] == he[base] );
			
			hePrev->next = heIndex+2;
			he[base]->next = heIndex+1;
#endif

			// Update the face indices of the half-edges to indicate they are in a triangle now.
			he[base]->face = faceIndex;

			m_vertData[f[prev]] = heIndex;
			m_faceData.push_back( heIndex );
			m_heData.push_back( newHe[0] );
			m_heData.push_back( newHe[1] );
			m_heData.push_back( newHe[2] );
			m_heData.push_back( newHe[3] );
		}
	}

	return faceIndex;
}

void EditMesh::delete_face( std::size_t f ){
	assert( f < m_faceData.size() );

	// We can assume that this face has 3 half-edges.
	std::size_t heIndices[3];
	heIndices[0] = m_faceData[f];
	
	half_edge* he[3];
	he[0] = &m_heData[heIndices[0]];
	he[1] = &m_heData[he[0]->next];
	he[2] = &m_heData[he[1]->next];

	heIndices[1] = he[0]->next;
	heIndices[2] = he[1]->next;
	
	assert( he[0]->face == f && he[1]->face == f && he[2]->face == f );
	assert( he[2]->next == m_faceData[f] );

	// Search for an edge that has a neighbor, but its prev edge doesn't. This is a canonical place to construct the algorithm from.
	std::size_t base = HOLE_INDEX;
	for( std::size_t i = 0; i < 3 && base == HOLE_INDEX; ++i ){
		if( m_heData[he[i]->twin].face != HOLE_INDEX && m_heData[he[(i+2)%3]->twin].face == HOLE_INDEX )
			base = i;
	}

	if( base == HOLE_INDEX ){
		if( m_heData[he[0]->twin].face == HOLE_INDEX ){
			// This is a lone triangle, so delete its half-edges and the exterior hole surrounding it too.

			// TODO: Remove the floating vertices? Currently we are leaving them.
			m_vertData[he[0]->vert] = HOLE_INDEX;
			m_vertData[he[1]->vert] = HOLE_INDEX;
			m_vertData[he[2]->vert] = HOLE_INDEX;

			// Delete all of the edges (both inside & outside half-edges). Must do this last since indices can change arbitrarily when deleting.
			std::size_t toDelete[] = { 
				heIndices[0], heIndices[1], heIndices[2], 
				he[0]->twin, he[1]->twin, he[2]->twin 
			};
			
			this->delete_half_edges_impl( toDelete );
			detail::delete_face( m_faceData, m_heData, f );
		}else{
			// This is an interior triangle. Only have to change the face_index to HOLE_INDEX for these edges.

			// Adjust any vertex references to new edges in non-hole faces.
			if( m_vertData[he[0]->vert] == heIndices[0] )
				m_vertData[he[0]->vert] = he[2]->twin;
			if( m_vertData[he[1]->vert] == heIndices[1] )
				m_vertData[he[1]->vert] = he[0]->twin;
			if( m_vertData[he[2]->vert] == heIndices[2] )
				m_vertData[he[2]->vert] = he[1]->twin;

			// Flag all these half-edges as being a hole now.
			he[0]->face = he[1]->face = he[2]->face = HOLE_INDEX;
			detail::delete_face( m_faceData, m_heData, f );
		}
	}else{
		std::rotate( he, he+base, he+3 );
		std::rotate( heIndices, heIndices+base, heIndices+3 );
		assert( m_heData[he[0]->twin].face != HOLE_INDEX );
		assert( m_heData[he[2]->twin].face == HOLE_INDEX );

		if( m_heData[he[1]->twin].face != HOLE_INDEX ){
			// We have one edge to remove, and a hole to connect to.
#ifdef USE_PREV
			he[1]->next = m_heData[he[2]->twin].next;
			he[0]->prev = m_heData[he[2]->twin].prev;
			m_heData[he[1]->next].prev = heIndices[1];
			m_heData[he[0]->prev].next = heIndices[0];
#else
			he[1]->next = m_heData[he[2]->twin].next;

			std::size_t hePrev = he[1]->next;
			while( m_heData[hePrev].next != he[2]->twin )
				hePrev = m_heData[hePrev].next;

			assert( m_heData[hePrev].next == he[2]->twin );
			m_heData[hePrev].next = heIndices[0];
#endif

			assert( m_heData[ m_vertData[ he[0]->vert ] ].face != HOLE_INDEX );
			assert( m_heData[ m_vertData[ he[1]->vert ] ].face != HOLE_INDEX );
			assert( m_heData[ m_vertData[ he[2]->vert ] ].face != HOLE_INDEX );

			// We may need to update the vertices if they referenced the edges we are deleting. Choose new half-edges that are inside 
			// non-hole triangles.
			if( m_vertData[he[0]->vert] == heIndices[0] )
				m_vertData[he[0]->vert] = m_heData[he[0]->twin].next;
			if( m_vertData[he[1]->vert] == heIndices[1] )
				m_vertData[he[1]->vert] = he[0]->twin;
			if( m_vertData[he[2]->vert] == heIndices[2] )
				m_vertData[he[2]->vert] = he[1]->twin;

			assert( m_heData[ m_vertData[ he[0]->vert ] ].face != HOLE_INDEX );
			assert( m_heData[ m_vertData[ he[1]->vert ] ].face != HOLE_INDEX );
			assert( m_heData[ m_vertData[ he[2]->vert ] ].face != HOLE_INDEX );

			he[0]->face = he[1]->face = HOLE_INDEX;

			std::size_t toDelete[] = { heIndices[2], he[2]->twin };

			// Delete the edges and face. Must do this last since indices can change arbitrarily when deleting.
			this->delete_half_edges_impl( toDelete );
			detail::delete_face( m_faceData, m_heData, f );
		}else{
			// We have two edges to remove, a vertex that will become floating, and a hole to connect to.
#ifdef USE_PREV
			he[0]->next = m_heData[he[1]->twin].next;
			he[0]->prev = m_heData[he[2]->twin].prev;
			m_heData[he[0]->next].prev = heIndices[0];
			m_heData[he[0]->prev].next = heIndices[0];
#else
			he[0]->next = m_heData[he[1]->twin].next;

			std::size_t hePrev = he[0]->next;
			while( m_heData[hePrev].next != he[2]->twin )
				hePrev = m_heData[hePrev].next;

			assert( m_heData[hePrev].next == he[2]->twin );
			m_heData[hePrev].next = heIndices[0];
#endif

			// We may need to update the vertices if they referenced the edges we are deleting. Choose new half-edges that are inside 
			// non-hole triangles.
			if( m_vertData[he[1]->vert] == heIndices[1] )
				m_vertData[he[1]->vert] = he[0]->twin;
			if( m_vertData[he[0]->vert] == he[2]->twin || m_vertData[he[0]->vert] == heIndices[0] )
				m_vertData[he[0]->vert] = m_heData[he[0]->twin].next;
			m_vertData[he[2]->vert] = HOLE_INDEX;

			// Update the face association of the one half-edge we are keeping (it joins the hole).
			he[0]->face = HOLE_INDEX;
			
			// Delete the edges and face. Must do this last since indices can change arbitrarily when deleting.
			std::size_t toDelete[] = { 
				heIndices[1], heIndices[2], 
				he[1]->twin, he[2]->twin 
			};

			this->delete_half_edges_impl( toDelete );
			detail::delete_face( m_faceData, m_heData, f );
		}
	}
}

std::size_t EditMesh::split_face_center( std::size_t f, std::size_t (*pOutFaceIndices)[3] ){
	assert( f < m_faceData.size() );
	assert( m_faceData[f] < m_heData.size() && m_heData[m_faceData[f]].face == f );

	std::size_t he[3];
	he[0] = m_faceData[f];
	he[1] = m_heData[he[0]].next;
	he[2] = m_heData[he[1]].next;

	assert( m_heData[he[2]].next == he[0] );
	assert( m_heData[he[0]].vert < m_vertices.size() && m_heData[he[1]].vert < m_vertices.size() && m_heData[he[2]].vert < m_vertices.size() );
	
	// New vert at face center
	Eigen::Vector3d newVert = ( m_vertices[ m_heData[ he[0] ].vert ] + m_vertices[ m_heData[ he[1] ].vert ] + m_vertices[ m_heData[ he[2] ].vert ] ) / 3.0;

	std::size_t newVertIndex = m_vertices.size();
	m_vertices.push_back( newVert );

	// Each half-edge gets associated to a new face, and we add 6 half-edges from the old vertices to the new.
	std::size_t newHeIndex = m_heData.size();
	std::size_t newFaceIndex = m_faceData.size();

	if( pOutFaceIndices ){
		(*pOutFaceIndices)[0] = f;
		(*pOutFaceIndices)[1] = newFaceIndex;
		(*pOutFaceIndices)[2] = newFaceIndex+1;
	}

	// Create six new half-edges connecting the center vertex to the old triangle corners.
	half_edge newHe[6];
	detail::init( newHe[0], newHeIndex+1, newHeIndex+3, m_heData[he[1]].vert, f );
	detail::init( newHe[1], he[0]       , newHeIndex+4, newVertIndex        , f );
	detail::init( newHe[2], newHeIndex+3, newHeIndex+5, m_heData[he[2]].vert, newFaceIndex );
	detail::init( newHe[3], he[1]       , newHeIndex  , newVertIndex        , newFaceIndex );
	detail::init( newHe[4], newHeIndex+5, newHeIndex+1, m_heData[he[0]].vert, newFaceIndex+1 );
	detail::init( newHe[5], he[2]       , newHeIndex+2, newVertIndex        , newFaceIndex+1 );

	// Connect the old half-edges to the new ones, and update their face association.
	//m_heData[he[0]].face = f;
	m_heData[he[0]].next = newHeIndex;
	m_heData[he[1]].face = newFaceIndex;
	m_heData[he[1]].next = newHeIndex+2;
	m_heData[he[2]].face = newFaceIndex+1;
	m_heData[he[2]].next = newHeIndex+4;

#ifdef USE_PREV
	newHe[0].prev = he[0];
	newHe[1].prev = newHeIndex;
	newHe[2].prev = he[1];
	newHe[3].prev = newHeIndex+2;
	newHe[4].prev = he[2];
	newHe[5].prev = newHeIndex+4;

	m_heData[he[0]].prev = newHeIndex+1;
	m_heData[he[1]].prev = newHeIndex+3;
	m_heData[he[2]].prev = newHeIndex+5;
#endif

	m_vertData.push_back( newHeIndex+3 ); // Arbitrary from 1, 3 & 5
	m_faceData[f] = he[0];
	m_faceData.push_back( he[1] );
	m_faceData.push_back( he[2] );
	m_heData.push_back( newHe[0] );
	m_heData.push_back( newHe[1] );
	m_heData.push_back( newHe[2] );
	m_heData.push_back( newHe[3] );
	m_heData.push_back( newHe[4] );
	m_heData.push_back( newHe[5] );

	return newVertIndex;
}

void EditMesh::split_boundary_edge( std::size_t heToSplit, std::size_t (*pOutVertIndices)[2], std::size_t (*pOutFaceIndices)[3] ){
	assert( heToSplit < m_heData.size() );
	assert( m_heData[heToSplit].face != HOLE_INDEX && m_heData[m_heData[heToSplit].twin].face == HOLE_INDEX );

	half_edge& heBase = m_heData[heToSplit];
	half_edge& heNext = m_heData[heBase.next];
	half_edge& hePrev = m_heData[heNext.next];
	half_edge& heTwin = m_heData[heBase.twin];

	Eigen::Vector3d newVert1 = m_vertices[ heBase.vert ] + ( m_vertices[ heNext.vert ] - m_vertices[ heBase.vert ] ) / 3.0;
	Eigen::Vector3d newVert2 = m_vertices[ heBase.vert ] + ( m_vertices[ heNext.vert ] - m_vertices[ heBase.vert ] ) * (2.0 / 3.0);

	// Construct 2 new faces and 8 new half_edges connecting the new verts to the off-edge vert.
	std::size_t newVertIndex = m_vertices.size();
	std::size_t newHeIndex = m_heData.size();
	std::size_t newFaceIndex = m_faceData.size();

	if( pOutVertIndices ){
		(*pOutVertIndices)[0] = newVertIndex;
		(*pOutVertIndices)[1] = newVertIndex+1;
	}

	if( pOutFaceIndices ){
		(*pOutFaceIndices)[0] = heBase.face;
		(*pOutFaceIndices)[1] = newFaceIndex;
		(*pOutFaceIndices)[2] = newFaceIndex+1;
	}

	half_edge newHe[8];
	detail::init( newHe[0], heNext.next , newHeIndex+1, newVertIndex   , heBase.face );
	detail::init( newHe[1], newHeIndex+4, newHeIndex  , hePrev.vert    , newFaceIndex );
	detail::init( newHe[2], newHeIndex+1, newHeIndex+3, newVertIndex+1 , newFaceIndex );
	detail::init( newHe[3], newHeIndex+5, newHeIndex+2, hePrev.vert    , newFaceIndex+1 );
	detail::init( newHe[4], newHeIndex+2, newHeIndex+6, newVertIndex   , newFaceIndex );
	detail::init( newHe[5], heBase.next , heBase.twin , newVertIndex+1 , newFaceIndex+1 );
	detail::init( newHe[6], newHeIndex+7, newHeIndex+4, newVertIndex+1 , HOLE_INDEX );
	detail::init( newHe[7], heTwin.next , heToSplit   , newVertIndex   , HOLE_INDEX );

#ifdef USE_PREV
	newHe[0].prev = heToSplit;
	newHe[1].prev = newHeIndex+2;
	newHe[2].prev = newHeIndex+4;
	newHe[3].prev = heBase.next;
	newHe[4].prev = newHeIndex+1;
	newHe[5].prev = newHeIndex+3;
	newHe[6].prev = heBase.twin;
	newHe[7].prev = newHeIndex+6;

	heNext.prev = newHeIndex+5;
	m_heData[heTwin.next].prev = newHeIndex+7
#endif

	heBase.next = newHeIndex;
	heBase.twin = newHeIndex+7;
	heTwin.next = newHeIndex+6;
	heTwin.twin = newHeIndex+5;
	heNext.next = newHeIndex+3;
	heNext.face = newFaceIndex+1;

	m_vertices.push_back( newVert1 );
	m_vertices.push_back( newVert2 );
	m_vertData.push_back( newHeIndex+4 );
	m_vertData.push_back( newHeIndex+5 );
	m_faceData[heBase.face] = heToSplit;
	m_faceData.push_back( newHeIndex+4 );
	m_faceData.push_back( newHeIndex+5 );
	m_heData.push_back( newHe[0] );
	m_heData.push_back( newHe[1] );
	m_heData.push_back( newHe[2] );
	m_heData.push_back( newHe[3] );
	m_heData.push_back( newHe[4] );
	m_heData.push_back( newHe[5] );
	m_heData.push_back( newHe[6] );
	m_heData.push_back( newHe[7] );
}

double EditMesh::get_cotan_weight( const vvert_iterator& it ) const {
	const half_edge *itCur = it.m_cur;
	const half_edge *itTwin = &m_heData[ itCur->twin ];

	double result = 0;

	Eigen::Vector3d a = this->get_vertex( itTwin->vert );
	Eigen::Vector3d b = this->get_vertex( itCur->vert );

	assert( ( itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX ) && "Invalid mesh: edge with no face on either side" );

	if( itCur->face != HOLE_INDEX ){
		Eigen::Vector3d c = this->get_vertex( this->prev( *itCur ).vert );
		Eigen::Vector3d e0 = b - c;
		Eigen::Vector3d e1 = a - c;

		// We use the dot product and norm of the cross product to get cos and sin respectively. cotan = cos / sin
		result += static_cast<double>( e0.dot(e1) ) / e0.cross(e1).norm();
	}

	if( itTwin->face != HOLE_INDEX ){
		Eigen::Vector3d c = this->get_vertex( this->prev( *itTwin ).vert );
		Eigen::Vector3d e0 = a - c;
		Eigen::Vector3d e1 = b - c;

		result += static_cast<double>( e0.dot(e1) ) / e0.cross(e1).norm();
	}
	
	return result;
}



Eigen::Vector3d EditMesh::get_normal( const vface_iterator& it ) const {
	Eigen::Vector3d a = this->get_vertex( it.m_next->vert );
	Eigen::Vector3d b = this->get_vertex( it.m_cur->vert );
	Eigen::Vector3d c = this->get_vertex( m_heData[ it.m_cur->next ].vert );
	
	return ( a - c ).cross( b - c ).normalized();
}

Eigen::Vector4d EditMesh::get_plane( const vface_iterator& it ) const {
	Eigen::Vector3d a = this->get_vertex( it.m_next->vert );
	Eigen::Vector3d b = this->get_vertex( it.m_cur->vert );
	Eigen::Vector3d c = this->get_vertex( m_heData[ it.m_cur->next ].vert );
	Eigen::Vector3d n = ( a - c ).cross( b - c ).normalized();

	// Plane equation Ax + By + Cz + D = 0 -> n.dot( [x,y,z] ) - n.dot( c ) = 0
	return Eigen::Vector4d( n.x(), n.y(), n.z(), -n.dot( c ) );
}

Eigen::Matrix4d get_quadric_error_matrix( EditMesh& m, std::size_t vertex ){
	Eigen::Matrix4d result = Eigen::Matrix4d::Zero();
	
	vface_iterator it;
	if( m.init_iterator( it, vertex ) ){
		do{
			// Skip holes. There can be only one per 1-ring or else this mesh is non-manifold.
			if( m.deref_iterator( it ) != HOLE_INDEX ){
				Eigen::Vector4d p = m.get_plane( it );
				result = ( result + p * p.transpose() ); // Construct the 4x4 Rank-1 distance matrix via outer product and add it in.
			}
		}while( m.advance_iterator( it ) );
	}

	return result;
}

double get_quadric_error( EditMesh& m, std::size_t v1, std::size_t v2 ){
	Eigen::Matrix4d errorMatrix = get_quadric_error_matrix( m, v1 ) + get_quadric_error_matrix( m, v2 );
	Eigen::Matrix4d derivMatrix = errorMatrix;
	derivMatrix.row( 3 ) = Eigen::Vector4d( 0, 0, 0, 1 );

	Eigen::FullPivLU< Eigen::Matrix4d > lu( derivMatrix );
	if( lu.isInvertible() ){
		// Find the point which minimizes the square error.
		Eigen::Vector4d	hp = lu.solve( Eigen::Vector4d(0,0,0,1) );
		return hp.dot( errorMatrix * hp );
	}else{
		// Use the edge midpoint.
		Eigen::Vector3d p = 0.5 * ( m.get_vertex( v1 ) + m.get_vertex( v2 ) );
		Eigen::Vector4d	hp( p.x(), p.y(), p.z(), 1.0 );
		return hp.dot( errorMatrix * hp );
	}
}

void collapse_one_edge( EditMesh& m )
{
	if( m.m_simplifyData.size() != m.m_vertData.size() )
	{
		for( std::size_t i = m.m_simplifyData.size(), iEnd = m.get_vert_size(); i < iEnd; ++i )
			m.m_simplifyData.push_back( get_quadric_error_matrix( m, i ) );
		
		m.m_simplifyQueue.resize( m.m_heData.size() );

		for( std::size_t i = 0, iEnd = m.m_heData.size(); i < iEnd; ++i )
		{
			// Only store data for one of the half_edges.
			if( i > m.m_heData[i].twin )
			{
				m.m_simplifyQueue[i].value = m.m_simplifyQueue[ m.m_heData[i].twin ].value;
			}
			else
			{
				std::size_t v1 = m.m_heData[ i ].vert;
				std::size_t v2 = m.m_heData[ m.m_heData[i].twin ].vert;

				m.m_simplifyQueue[ i ].value = get_quadric_error( m, v1, v2 );
			}
		}
	}

	Eigen::Matrix4d errorMatrix;
	std::size_t v;

	do{
		double bestValue = std::numeric_limits<double>::max();
		std::size_t bestIndex = HOLE_INDEX;

		// Find the lowest error half_edge. This *could* be optimized with a priority queue or heap, but that's freaking complicated ...
		for( std::size_t i = 0, iEnd = m.m_heData.size(); i < iEnd; ++i )
		{
			if( m.m_simplifyQueue[i].value < bestValue )
			{
				bestValue = m.m_simplifyQueue[i].value;
				bestIndex = i;
			}
		}

		if( bestIndex == HOLE_INDEX )
		{
			// We couldn't collapse anything.
			std::cerr << "Can't collapse anything!" << std::endl;
			return;
		}

		// Store the new error matrix while our vert indices are still valid.
		errorMatrix = m.m_simplifyData[ m.m_heData[bestIndex].vert ] + m.m_simplifyData[ m.m_heData[ m.m_heData[bestIndex].twin ].vert ];

		v = m.collapse_edge( bestIndex );

		if( v == HOLE_INDEX )
			m.m_simplifyQueue[ bestIndex ].value = std::numeric_limits<double>::quiet_NaN();

	 }while( v == HOLE_INDEX );

	m.m_simplifyData[v] = errorMatrix;

	// Update the error values for all the verts in the new 1-ring.
	vvert_iterator it;
	if( m.init_iterator( it, v ) ){
		do{
			std::size_t he = m.deref_iterator_left_edge( it );
			std::size_t heTwin = m.deref_iterator_right_edge( it );

			m.m_simplifyQueue[ he ].value = m.m_simplifyQueue[ heTwin ].value = get_quadric_error( m, v, m.deref_iterator( it ) );
		}while( m.advance_iterator( it ) );
	}

	m.edit_count_updater();
}

Eigen::Vector3d EditMesh::get_vnormal( std::size_t i ) const {
    vvert_iterator vit;
    init_iterator(vit, i);

    Eigen::Vector3d normal;
    normal.setZero();

    Eigen::Vector3d center = this->get_vertex( i );
    Eigen::Vector3d vec_prev;
    Eigen::Vector3d vec_curr = this->get_vertex( deref_iterator(vit) ) - center;

    advance_iterator(vit);
    vit.m_end = vit.m_cur;
    do {
        vec_prev = vec_curr;
        vec_curr = this->get_vertex( deref_iterator(vit) ) - center;

        if (m_heData[vit.m_cur->twin].face != HOLE_INDEX)
            normal += vec_curr.cross(vec_prev);
    } while (advance_iterator(vit));

    return normal.normalized();
}

Eigen::Vector3d EditMesh::get_fnormal( std::size_t i ) const {
	const half_edge *e1 = &m_heData[ m_faceData[ i ] ];
    const half_edge *e2 = &this->next(*e1);
    const half_edge *e3 = &this->next(*e2);

    Eigen::Vector3d a = this->get_vertex( e1->vert );
    Eigen::Vector3d b = this->get_vertex( e2->vert );
    Eigen::Vector3d c = this->get_vertex( e3->vert );

    Eigen::Vector3d normal = (b - a).cross(c - a);

    return normal.normalized();
}

void EditMesh::getIndicesForFace(size_t tri_index, size_t indicesForFace[3]) {
	fvert_iterator fvit;
	init_iterator( fvit, tri_index );

	for( size_t i = 0; i < 3; i++ ) {
		indicesForFace[i] = deref_iterator(fvit);
		advance_iterator( fvit );
	}
}

Eigen::Vector3d EditMesh::getFaceMidpoint(size_t tri_index) {
	fvert_iterator fvit;
	init_iterator( fvit, tri_index );

    Eigen::Vector3d average = Eigen::Vector3d::Zero();
    int count = 0;
	for( size_t i = 0; i < 3; i++ ) {
        size_t index = deref_iterator(fvit);
		if (index != HOLE_INDEX) {
            count++;
            average += get_vertex(index);
        }
		advance_iterator( fvit );
	}
    return average / count;
}

void EditMesh::get_draw_data( float *verts, int *indices ) const {

    /* get each vertex only once. This is good for efficiency
     * but results in bad looking meshes due to each vertex 
     * having a fixed normal

        for( std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; i++ ){
            const half_edge* he = &m_heData[ m_faceData[i] ];

            for( int j = 0; j < 3; j++){
                indices[3*i +j] = he->vert;
                he = &this->next(*he);
            }
        }

        for( std::size_t i = 0, iEnd = m_vertData.size(); i < iEnd; i++ ){
		    Eigen::Vector3d vert = this->get_vertex( i );
            for( int j = 0; j < 3; j++)
                verts[3*i+j] = (float) vert[j];
        }
    */

    // for each face
    for( std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; i++ ){
        const half_edge* he = &m_heData[ m_faceData[i] ];

        // for each vertex of the face
        for( int j = 0; j < 3; j++){
            Eigen::Vector3d vert = this->get_vertex(he->vert);
            indices[3*i+j] = 3*i+j;

            // for each component of the vertex
            for( int k = 0; k < 3; k++){
                verts[3*(3*i+j) + k] = vert[k];
            }
            he = &this->next(*he);
        }
    }
}

void EditMesh::get_draw_normals( float *normals ) const {

    /* this finds the averaged vertex normals which results in
     * poor looking meshes when they are not smooth

        for( std::size_t i = 0, iEnd = m_vertData.size(); i < iEnd; i++ ){
		    Eigen::Vector3d normal = this->get_normal( i );
            for( int j = 0; j < 3; j++)
                normals[3*i+j] = (float) normal[j];
        }
    */

    for( std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; i++ ){
		Eigen::Vector3d normal = this->get_fnormal( i );

        for( int j = 0; j < 3; j++){
            for( int k = 0; k < 3; k++){
                normals[3*(3*i+j) + k] = normal[k];
            }
        }
    }
}

// call instead of init to test edge flip
// easiest way is hacking it into mesh constructor
void EditMesh::test_flip() {
    std::vector<double> xyz;
    std::vector<std::size_t> faces;

    // four verts
    xyz.push_back(-1); xyz.push_back(0); xyz.push_back(0);
    xyz.push_back(0); xyz.push_back(1); xyz.push_back(1);
    xyz.push_back(0); xyz.push_back(1); xyz.push_back(-1);
    xyz.push_back(1); xyz.push_back(0); xyz.push_back(0);

    // two triangles
    faces.push_back(0); faces.push_back(1);
    faces.push_back(2); faces.push_back(2);
    faces.push_back(1); faces.push_back(3);

    this->init(xyz, faces);

    half_edge *he = &m_heData[0];
    for( int i = 0; i < m_heData.size(); i++ ) {
        he = &m_heData[i];
        if (he->face != HOLE_INDEX &&
            m_heData[ he->twin ].face != HOLE_INDEX )
            break;
    }
    flip_edge(*he);
    this->edit_count++;
}

void EditMesh::write_to_obj_stream( std::ostream& stream ) const {
	for( auto& v : m_vertices )
		stream << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << std::endl;
	stream << std::endl;
	for( std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; ++i ){
		fvert_iterator it;
		this->init_iterator( it, i );
		stream << "f ";
		bool isFirst = true;
		do{
			if( !isFirst )
				stream << ' ';
			isFirst = false;
			stream << this->deref_iterator( it )+1;
		}while( this->advance_iterator( it ) );
		stream << std::endl;
	}
}

void EditMesh::verify() const {
	for( std::size_t i = 0, iEnd = m_faceData.size(); i < iEnd; ++i ){
		std::size_t c = 0;
		
		const half_edge* it = &m_heData[ m_faceData[i] ];
		assert( it->next != m_faceData[i] );
		while( it->next != m_faceData[i] ){
			assert( it->face == i );
			assert( it->next != HOLE_INDEX && it->twin != HOLE_INDEX && it->vert < m_vertData.size() );
			assert( ( m_heData[ it->twin ].face == HOLE_INDEX || m_heData[ m_heData[it->next].twin ].face != m_heData[ it->twin ].face ) && "Can't have two edges shared between the same faces!" );
			it = &m_heData[it->next];
			assert( ++c < 1000000 ); // This isn't strictly a problem, but probably no face has a million verts in it.
		}
	}

	for( std::size_t i = 0, iEnd = m_vertData.size(); i < iEnd; ++i ){
		assert( m_vertData[i] == HOLE_INDEX || m_vertData[i] < m_heData.size() );
		if( m_vertData[i] != HOLE_INDEX ){
			const half_edge* it = &m_heData[ m_vertData[i] ];
			assert( it->vert == i );
			assert( it->face != HOLE_INDEX && "By convention, vertices should not reference hole faces" );
		}
	}

	for( std::size_t i = 0, iEnd = m_heData.size(); i < iEnd; ++i ){
		const half_edge* it = &m_heData[i];
		assert( it->vert < m_vertData.size() );
		assert( it->face == HOLE_INDEX || it->face < m_faceData.size() );

		assert( it->next < m_heData.size() );
		assert( it->next != i );
		assert( m_heData[it->next].face == it->face );
		assert( m_heData[it->next].vert != it->vert );

		assert( it->twin < m_heData.size() );
		assert( m_heData[it->twin].twin == i );
		assert( m_heData[it->twin].face != it->face );
		assert( m_heData[it->twin].vert == m_heData[it->next].vert );

#ifdef USE_PREV
		assert( it->prev < m_heData.size() );
		assert( it->prev != i );
		assert( m_heData[it->next].prev == i );
		assert( m_heData[it->prev].next == i );
		assert( m_heData[it->prev].face == it->face );
		assert( m_heData[it->prev].vert != it->vert );
#endif
	}
}

void EditMesh::test(){
	EditMesh m1, m2, m3;

	std::vector<double> v;
	std::vector<std::size_t> f;

	v.push_back( 0 ); v.push_back( 0 ); v.push_back( 0 );
	v.push_back( 1 ); v.push_back( 0 ); v.push_back( 0 );
	v.push_back( 0 ); v.push_back( 1 ); v.push_back( 0 );
	v.push_back( 1 ); v.push_back( 1 ); v.push_back( 0 );
	v.push_back( 2 ); v.push_back( 1 ); v.push_back( 0 );
	v.push_back( 1 ); v.push_back( 2 ); v.push_back( 0 );

	f.push_back( 0 ); f.push_back( 1 ); f.push_back( 2 );
	f.push_back( 2 ); f.push_back( 1 ); f.push_back( 3 );
	f.push_back( 3 ); f.push_back( 1 ); f.push_back( 4 );
	f.push_back( 4 ); f.push_back( 5 ); f.push_back( 3 );
	f.push_back( 3 ); f.push_back( 5 ); f.push_back( 2 );

	m1.init( v, f );

	for( std::size_t i = 0, iEnd = v.size(); i < iEnd; i += 3 )
		assert( m2.add_vertex( v[i], v[i+1], v[i+2] ) == i/3 );
	for( std::size_t i = 0, iEnd = f.size(); i < iEnd; i += 3 )
		assert( m2.add_face( *reinterpret_cast<std::size_t(*)[3]>( &f[i] ) ) == i/3 );

	assert( m1.get_face_size() == m2.get_face_size() );
	assert( m1.get_vert_size() == m2.get_vert_size() );

	m1.verify();
	m2.verify();

	m2.delete_face( 0 );
	m2.verify();

	m2.delete_face( 0 ); // Was face 4 originally
	m2.verify();

	m2.delete_face( 0 ); // Was face 3 originally
	m2.verify();

	m2.delete_face( 0 ); // Was face 2 originally
	m2.verify();

	m2.delete_face( 0 ); // Was face 1 originally
	m2.verify();

	assert( m2.get_face_size() == 0 );

	m2.add_face( 0, 1, 2 );
	m2.verify();

	m2.split_face_center( 0 );
	m2.verify();

	assert( m2.get_face_size() == 3 );

	m2.clear();
	m2.add_vertex( 0, 0, 0 );
	m2.add_vertex( 1, 0, 0 );
	m2.add_vertex( 0, 1, 0 );
	m2.add_face( 0, 1, 2 );
	m2.split_boundary_edge( m2.find_twin( 0, 1 )->twin );
	m2.verify();

	assert( m2.get_face_size() == 3 );

	m3.add_vertex( -1, 0, 0 );
	m3.add_vertex( 1, 0, 0 );
	m3.add_vertex( 0, 1, 0 );
	m3.add_vertex( 0, -1, 0 );
	m3.add_vertex( -1, 1, 0 );
	m3.add_vertex( -1, -1, 0 );
	m3.add_vertex( 1, 1, 0 );
	m3.add_vertex( 1, -1, 0 );

	m3.add_face( 0, 1, 2 );
	m3.add_face( 0, 2, 4 );
	m3.add_face( 0, 4, 5 );
	m3.add_face( 0, 5, 3 );
	m3.add_face( 0, 3, 1 );
	
	m3.add_face( 1, 3, 7 );
	m3.add_face( 1, 7, 6 );
	m3.add_face( 1, 6, 2 );
	m3.verify();

	m3.flip_edge( *m3.find_edge( 1, 0 ) );
	m3.verify();
	m3.flip_edge( *m3.find_edge( 2, 3 ) );
	m3.verify();

	std::size_t newVert = m3.collapse_edge( m3.find_edge( 1, 0 )->twin );
	assert( newVert != HOLE_INDEX );
	m3.verify();

	/*m2.clear();
	m2.add_vertex( 0.010744, 0.483695, 0.298761 );
	m2.add_vertex( 0.010538, 0.484281, 0.305409 );
	m2.add_vertex( 0.014906, 0.48369, 0.304997 );
	m2.add_vertex( 0.006473, 0.484811, 0.30548 );
	m2.add_vertex( 0.010333, 0.484867, 0.312038 );
	m2.add_vertex( 0.004998, 0.485704, 0.314376 );
	m2.add_vertex( 0.010129, 0.485783, 0.323883 );
	m2.add_vertex( 0.016209, 0.484307, 0.313866 );
	m2.add_face( 7, 1, 4 );
	m2.add_face( 7, 6, 1 );
	m2.add_face( 6, 3, 1 );
	m2.add_face( 5, 3, 6 );
	m2.add_face( 5, 0, 3 );
	m2.add_face( 3, 2, 4 );
	m2.add_face( 3, 0, 2 );
	m2.add_face( 1, 3, 4 );

	for( std::size_t i = m2.get_face_size(); i > 0; --i ){
		m2.delete_face( i-1 );
		m2.verify();
	}*/

	m1.clear();
	m1.add_vertex( -1, 0, 0 );
	m1.add_vertex( 1, 0, 0 );
	m1.add_vertex( 0, 1, 0 );
	m1.add_vertex( -2, 1, 0 );
	m1.add_vertex( -2, -1, 0 );
	m1.add_vertex( 0, -1, 0 );
	m1.add_vertex( 2, -1, 0 );
	m1.add_vertex( 2, 1, 0 );
	m1.add_vertex( -2, 2, 0 );
	m1.add_vertex( -3, 1, 0 );
	m1.add_vertex( -3, -1, 0 );
	m1.add_vertex( -2, -2, 0 );
	m1.add_vertex( 2, -2, 0 );
	m1.add_vertex( 3, -1, 0 );
	m1.add_vertex( 3, 1, 0 );
	m1.add_vertex( 2, 2, 0 );

	m1.add_face( 0, 1, 2 );
	m1.add_face( 0, 2, 3 );
	m1.add_face( 0, 3, 4 );
	m1.add_face( 0, 4, 5 );
	m1.add_face( 0, 5, 1 );
	m1.add_face( 1, 5, 6 );
	m1.add_face( 1, 6, 7 );
	m1.add_face( 1, 7, 2 );
	m1.add_face( 3, 2, 8 );
	m1.add_face( 3, 8, 9 );
	m1.add_face( 3, 9, 10 );
	m1.add_face( 3, 10, 4 );
	m1.add_face( 4, 10, 11 );
	m1.add_face( 4, 11, 5 );
	m1.add_face( 5, 11, 12 );
	m1.add_face( 5, 12, 6 );
	m1.add_face( 6, 12, 13 );
	m1.add_face( 6, 13, 14 );
	m1.add_face( 6, 14, 7 );
	m1.add_face( 7, 14, 15 );
	m1.add_face( 7, 15, 2 );
	m1.add_face( 2, 15, 8 );

	m1.verify();

	std::vector< std::size_t > faces;
	vface_iterator it;
	if( m1.init_iterator( it, 0 ) ){
		do{
			std::vector< std::size_t >::iterator itInsert = std::lower_bound( faces.begin(), faces.end(), m1.deref_iterator( it ), std::greater<std::size_t>() );
			if( itInsert == faces.end() || *itInsert != m1.deref_iterator( it ) )
				faces.insert( itInsert, m1.deref_iterator( it ) );
		}while( m1.advance_iterator( it ) );
	}
	if( m1.init_iterator( it, 1 ) ){
		do{
			std::vector< std::size_t >::iterator itInsert = std::lower_bound( faces.begin(), faces.end(), m1.deref_iterator( it ), std::greater<std::size_t>() );
			if( itInsert == faces.end() || *itInsert != m1.deref_iterator( it ) )
				faces.insert( itInsert, m1.deref_iterator( it ) );
		}while( m1.advance_iterator( it ) );
	}

	std::swap( faces[faces.size()-1], faces[faces.size()-2] );

	for( auto face : faces ){
		m1.delete_face( face );
		m1.verify();
	}

	std::unique_ptr<EditMesh> em( loadEditMeshFromFile("Mesh/Collapse-2Ring.obj") );
	em->verify();
	em->collapse_edge( 5, 4 );
	em->verify();

	/*em.reset( loadEditMeshFromFile("Mesh/camel.obj") );

	for( std::size_t i = 0; i < 9667; ++i ){
		collapse_one_edge( *em );
		if ( i > 9000 || i % 1000 == 0 )
			em->verify();
	}*/

	em->clear();
	em->add_vertex(2, .345, 0);
	em->add_vertex(4, 0, 0);
	em->add_vertex(-4, 0, 0);
	em->add_vertex(0, 1, 0);
	em->add_vertex(0, -1, 0);

	em->add_face(0, 3, 2);
	em->add_face(0, 2, 4);
	em->add_face(0, 4, 1);
	em->add_face(0, 1, 3);
}

bool EditMesh::is_safe_addface( std::size_t v1, std::size_t v2, std::size_t v3 ) {

    std::set<std::size_t> vv;
    vv.insert( m_vertData[v1] );
    vv.insert( m_vertData[v2] );
    vv.insert( m_vertData[v3] );


    // if there's one disconnected vertex its safe
    // more than one is not
    if( vv.size() < 3 )
        return false;
    else if( vv.count(HOLE_INDEX) > 0 )
        return true;

    vv.clear();
    vv.insert(v1);
    vv.insert(v2);
    vv.insert(v3);

    int count = 0;
    vvert_iterator vit;
    init_iterator(vit, v1);
    do {
        if (vit.m_cur->vert == v2 ||
            vit.m_cur->vert == v3) {
            count++;
            break;
        }
    } while( this->advance_iterator(vit) );

    init_iterator(vit, v2);
    do {
        if (vit.m_cur->vert == v1 ||
            vit.m_cur->vert == v3) {
            count++;
            break;
        }
    } while( this->advance_iterator(vit) );

    // if two of the vertices have a next within the triplet, its safe
    if ( count > 1 )
        return true;

    // else unsafe
    return false;
    //return (v1 != HOLE_INDEX && v2 != HOLE_INDEX) ||
    //       (v2 != HOLE_INDEX && v3 != HOLE_INDEX) ||
    //       (v3 != HOLE_INDEX && v1 != HOLE_INDEX);
}

void EditMesh::Subdivision(std::size_t subdivMethod)
{
	// is open mesh
	std::size_t counterBoundaryHalfEdge(0);
	for ( std::size_t i = 0, iEnd = m_heData.size(); i < iEnd; ++i )
	{
		if (m_heData[i].face == HOLE_INDEX)
			++counterBoundaryHalfEdge;
	}

	// Save current information
	std::vector< half_edge > tempHeData(m_heData);			// Store current half-edges that make up the mesh.
	std::vector< std::size_t > tempFaceData(m_faceData);	// Store current face index data
	std::vector< std::size_t > tempVertData(m_vertData);	// Store current vertex index data
    std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > tempVertices(m_vertices);

	// resize each container
	m_heData.resize( (tempHeData.size() -  counterBoundaryHalfEdge) * 4 + counterBoundaryHalfEdge * 2); 
	m_faceData.resize( tempFaceData.size() * 4);
	m_vertData.resize( tempVertData.size() + tempHeData.size() / 2); // exclude the boundary case, but vitual boundary still can be divided into 2
	m_vertices.resize( tempVertices.size() + tempHeData.size() / 2);
	
	// insert new vertices
	std::size_t newVertPosition(tempVertices.size());
	std::vector< bool > isVisited( tempHeData.size(), false);
	for ( std::size_t i = 0, iEnd = tempHeData.size(); i < iEnd; ++i )
	{
		// is visited or polluted?
		if (isVisited[i])
			continue;
		
		// Calculate the new vertex's position
		Eigen::Vector3d topLeft, top, topRight, left, right, lowLeft, low, lowRight, newVertex;
		std::size_t indexTopLeft, indexTop, indexTopRight, indexLeft, indexRight, indexLowLeft, indexLow, indexLowRight;
		half_edge hf(tempHeData[i]), hfTwin(tempHeData[tempHeData[i].twin]);
		if (hf.face == HOLE_INDEX)	// skip the boundary half edge case, it is processed when meet its twin half edge
		{
			continue;
		}
		else if (hfTwin.face == HOLE_INDEX)	// the boundary half edge's twin need special processing: four point scheme
		{
			indexLeft     = tempHeData[i].vert;
			indexRight	  = tempHeData[tempHeData[i].twin].vert;
			std::size_t indexLeftNext(tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].next].twin].next].twin].vert);
			std::size_t indexRightNext(tempHeData[tempHeData[i].next].vert);

			Eigen::Vector3d mid[] = {1.0 /2 * (tempVertices[indexLeft] + tempVertices[indexRight]), 
										1.0 /2 * (tempVertices[indexLeftNext] + tempVertices[indexRightNext])};
			
			newVertex = (13 * mid[0] - mid[1]) / 12.0;
		}
		else if (!isBoundaryFace(tempHeData, tempFaceData, hf.face) && !isBoundaryFace(tempHeData, tempFaceData, hfTwin.face))	// the interior case that is shared by open and closed cases
		{
			indexTopLeft  = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].next].next].twin].next].next].vert;
			indexTop      = tempHeData[tempHeData[tempHeData[i].next].next].vert;
			indexTopRight = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].next].twin].next].next].vert;
			indexLeft     = tempHeData[i].vert;
			indexRight	  = tempHeData[tempHeData[i].twin].vert;
			indexLowLeft  = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].twin].next].next].vert;
			indexLow      = tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].vert;
			indexLowRight = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].twin].next].next].vert;
			topLeft  = m_vertices[indexTopLeft];
			top		 = m_vertices[indexTop];
			topRight = m_vertices[indexTopRight];
			left     = m_vertices[indexLeft];
			right    = m_vertices[indexRight];
			lowLeft  = m_vertices[indexLowLeft];
			low      = m_vertices[indexLow];
			lowRight = m_vertices[indexLowRight];
		}
		else if (isBoundaryFace(tempHeData, tempFaceData, hf.face) && !isBoundaryFace(tempHeData, tempFaceData, hfTwin.face)) // hf's face is boundary face
		{
			indexTop      = tempHeData[tempHeData[tempHeData[i].next].next].vert;
			indexLeft     = tempHeData[i].vert;
			indexRight	  = tempHeData[tempHeData[i].twin].vert;
			indexLowLeft  = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].twin].next].next].vert;
			indexLow      = tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].vert;
			indexLowRight = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].twin].next].next].vert;
			top		 = m_vertices[indexTop];
			left     = m_vertices[indexLeft];
			right    = m_vertices[indexRight];
			lowLeft  = m_vertices[indexLowLeft];
			low      = m_vertices[indexLow];
			lowRight = m_vertices[indexLowRight];

			// get the topLeft and topRight by refelecting
			topLeft  = left + top - right;
			topRight = top + right - left;
		}
		else if (isBoundaryFace(tempHeData, tempFaceData, hf.face) && isBoundaryFace(tempHeData, tempFaceData, hfTwin.face))
		{
			indexTop      = tempHeData[tempHeData[tempHeData[i].next].next].vert;
			indexLeft     = tempHeData[i].vert;
			indexRight	  = tempHeData[tempHeData[i].twin].vert;
			indexLow      = tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].vert;
			top		 = m_vertices[indexTop];
			left     = m_vertices[indexLeft];
			right    = m_vertices[indexRight];
			low      = m_vertices[indexLow];

			// get the topLeft and topRight by reflecting
			topLeft  = left + top - right;
			topRight = top + right - left;
			lowLeft  = left +  low - right;
			lowRight = low  + right - left;
		}
		else if (!isBoundaryFace(tempHeData, tempFaceData, hf.face) && isBoundaryFace(tempHeData, tempFaceData, hfTwin.face))
		{
			indexTopLeft  = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].next].next].twin].next].next].vert;
			indexTop      = tempHeData[tempHeData[tempHeData[i].next].next].vert;
			indexTopRight = tempHeData[tempHeData[tempHeData[tempHeData[tempHeData[i].next].twin].next].next].vert;
			indexLeft     = tempHeData[i].vert;
			indexRight	  = tempHeData[tempHeData[i].twin].vert;
			indexLow      = tempHeData[tempHeData[tempHeData[tempHeData[i].twin].next].next].vert;
			topLeft  = m_vertices[indexTopLeft];
			top		 = m_vertices[indexTop];
			topRight = m_vertices[indexTopRight];
			left     = m_vertices[indexLeft];
			right    = m_vertices[indexRight];
			low      = m_vertices[indexLow];

			// get the topLeft and topRight by reflecting
			lowLeft  = left +  low - right;
			lowRight = low  + right - left;
		}

		if (hfTwin.face != HOLE_INDEX)
		{
			if (subdivMethod)	
				newVertex = 1.0 / 8 * ( top + low ) + 3.0 / 8 * ( left + right );
			else newVertex = -1.0 / 16 * ( topLeft + topRight + lowLeft + lowRight ) + 2.0 / 16 * ( top + low ) + 8.0 / 16 * ( left + right );
		}
		m_vertices[newVertPosition] = newVertex;
		
		// new half edges
		std::size_t index[] = {2 * i, 2 * i + 1, 2 * tempHeData[i].twin, 2 * tempHeData[i].twin + 1};
		detail::init( m_heData[index[0]], HOLE_INDEX, index[3], indexLeft, HOLE_INDEX );
		detail::init( m_heData[index[1]], HOLE_INDEX, index[2], newVertPosition, HOLE_INDEX );
		detail::init( m_heData[index[2]], HOLE_INDEX, index[1], indexRight, HOLE_INDEX );
		detail::init( m_heData[index[3]], HOLE_INDEX, index[0], newVertPosition, HOLE_INDEX );

		// assign the new and old vertices with two half edges respectively
		m_vertData[indexLeft] = index[0];
		m_vertData[indexRight] = index[2];
		m_vertData[newVertPosition] = index[1];
		// end of assigning

		++newVertPosition;
		isVisited[i] = isVisited[tempHeData[i].twin] = true;	// The twin has already been visited.
	}
	// end of inserting new vertices

	for ( std::size_t i = 0, iEnd = isVisited.size(); i < iEnd; ++i )
	{
		if (!isVisited[i])
			std::cerr << "error occurred when insert new vertices." << std::endl;

		// for reuse of this flag vector
		isVisited[i] = false;
	}


	// Mesh structure update
	std::size_t newHEdgePosition(tempHeData.size() * 2);
	for ( std::size_t i = 0, iEnd = tempHeData.size(); i < iEnd; ++i )
	{
		// is visited or polluted?
		if ( isVisited[i] )
			continue;

		// Complete half edges for the triangle the half edge m_heData[i] lies.
		std::size_t heIndexOld[] = {i, tempHeData[i].next, tempHeData[tempHeData[i].next].next};
		std::size_t heIndexNew[] = {2 * heIndexOld[0], 2 * heIndexOld[0] + 1, newHEdgePosition, newHEdgePosition + 1, 
										2 * heIndexOld[1], 2 * heIndexOld[1] + 1, newHEdgePosition + 2, newHEdgePosition + 3,
											2 * heIndexOld[2], 2 * heIndexOld[2] + 1, newHEdgePosition + 4, newHEdgePosition + 5};

		if (tempHeData[i].face != HOLE_INDEX)	// for boundary half edges, there is no four half edges, only two
		{
			// for the first edge
			detail::init( m_heData[heIndexNew[2]], HOLE_INDEX, heIndexNew[3], m_heData[heIndexNew[9]].vert, HOLE_INDEX );
			detail::init( m_heData[heIndexNew[3]], HOLE_INDEX, heIndexNew[2], m_heData[heIndexNew[5]].vert, HOLE_INDEX );
			// for the next edge
			detail::init( m_heData[heIndexNew[6]], HOLE_INDEX, heIndexNew[7], m_heData[heIndexNew[1]].vert, HOLE_INDEX );
			detail::init( m_heData[heIndexNew[7]], HOLE_INDEX, heIndexNew[6], m_heData[heIndexNew[9]].vert, HOLE_INDEX );
			// for the next next edge
			detail::init( m_heData[heIndexNew[10]], HOLE_INDEX, heIndexNew[11], m_heData[heIndexNew[5]].vert, HOLE_INDEX );
			detail::init( m_heData[heIndexNew[11]], HOLE_INDEX, heIndexNew[10], m_heData[heIndexNew[1]].vert, HOLE_INDEX );
		}
		// end of complete

		// order these half edges (i.e., next edge) in an anticlockwise way
		if (tempHeData[i].face != HOLE_INDEX)
		{			
			m_heData[heIndexNew[0]].next = heIndexNew[6];
			m_heData[heIndexNew[1]].next = heIndexNew[4];
			m_heData[heIndexNew[2]].next = heIndexNew[5];
			m_heData[heIndexNew[3]].next = heIndexNew[7];
			m_heData[heIndexNew[4]].next = heIndexNew[10];
			m_heData[heIndexNew[5]].next = heIndexNew[8];
			m_heData[heIndexNew[6]].next = heIndexNew[9];
			m_heData[heIndexNew[7]].next = heIndexNew[11];
			m_heData[heIndexNew[8]].next = heIndexNew[2];
			m_heData[heIndexNew[9]].next = heIndexNew[0];
			m_heData[heIndexNew[10]].next = heIndexNew[1];
			m_heData[heIndexNew[11]].next = heIndexNew[3];
		} 
		else
		{
			m_heData[heIndexNew[0]].next = heIndexNew[1];
			m_heData[heIndexNew[1]].next = heIndexNew[4];
		}
		// end of ordering

		if (tempHeData[i].face != HOLE_INDEX)
		{
			// add half edge map for m_faceData
			std::size_t faceIndex(tempHeData[i].face);
			m_faceData[4 * faceIndex + 0] = heIndexNew[0];
			m_faceData[4 * faceIndex + 1] = heIndexNew[4];
			m_faceData[4 * faceIndex + 2] = heIndexNew[8];
			m_faceData[4 * faceIndex + 3] = heIndexNew[3];	
			// end of adding
	
			// add face information for half edges
			m_heData[heIndexNew[0]].face = m_heData[heIndexNew[6]].face = m_heData[heIndexNew[9]].face = 4 * faceIndex + 0;
			m_heData[heIndexNew[1]].face = m_heData[heIndexNew[4]].face = m_heData[heIndexNew[10]].face = 4 * faceIndex + 1;
			m_heData[heIndexNew[5]].face = m_heData[heIndexNew[8]].face = m_heData[heIndexNew[2]].face = 4 * faceIndex + 2;
			m_heData[heIndexNew[3]].face = m_heData[heIndexNew[7]].face = m_heData[heIndexNew[11]].face = 4 * faceIndex + 3;
			// end of adding

			newHEdgePosition += 6;
			isVisited[heIndexOld[0]] = isVisited[heIndexOld[1]] = isVisited[heIndexOld[2]] = true;	// The face has already been visited.
		}

		isVisited[heIndexOld[0]] = true;
	}

	for ( std::size_t i = 0, iEnd = isVisited.size(); i < iEnd; ++i )
	{
		if (!isVisited[i])
			std::cerr << "error occurred when update mesh structure." << std::endl;

		// for reuse of this flag vector
		isVisited[i] = false;
	}

	// update the old vertices for loop scheme
	if (subdivMethod)
	{
		for ( std::size_t i = 0, iEnd = tempVertData.size(); i < iEnd; ++i )
		{
			vvert_iterator vvIter;
			if (!this->init_iterator( vvIter, i))
				std::cerr << "there is a float point when updating old vertices by loop scheme" << std::endl;
			
			if (this->reset_boundary_iterator(vvIter))		// boundary vertices need special processing
			{
				Eigen::Vector3d left(m_vertices[vvIter.m_cur->vert]);
				Eigen::Vector3d right(m_vertices[m_heData[m_heData[vvIter.m_cur->next].twin].vert]);

				m_vertices[i] = 1.0 / 2 * (left + right);

			}
			else
			{
				std::size_t valence(0);
				Eigen::Vector3d sum1Ring(0, 0, 0);
					
				// sum the 1 ring neighbors
				do 
				{
					sum1Ring += tempVertices[tempHeData[m_heData[vvIter.m_cur->twin].twin / 2].vert];
					++valence;
				} while (this->advance_iterator(vvIter));
				// end summing
			
				// weight for the center vertex
				double centerWeight(64.0 * valence / (40 - (3 + 2 * std::cos(2 * PI / valence)) * (3 + 2 * std::cos(2 * PI / valence))) - valence);
			
				// update center vertex
				m_vertices[i] = (sum1Ring + centerWeight * m_vertices[i]) / (valence + centerWeight);
			}

		}
	}


	this->edit_count_updater();	// for mesh view
}

void EditMesh::get_color_vert(vector<size_t>& _vertVecColor)
{
	_vertVecColor = vertVecColor;
}

/*========================================
* error computation functions 
*========================================*/
double EditMesh::get_Gauss_error(std::size_t vertex )
{
	double sum(0);

	// Boundary handling
	if (this->isBoundaryVert(vertex))
	{
		//std::cerr << "error because of boundary point" << std::endl;
		return std::numeric_limits<double>::max();
	}

	// initial the iterator visit the one ring of "vertex"
	vvert_iterator it;
	if (!this->init_iterator(it, vertex))
		std::cerr << "error because of float point" << std::endl;

	// sum the 1-ring angles	
	do {
		Eigen::Vector3d v1 = m_vertices[it.m_cur->vert] - m_vertices[vertex];
		half_edge he = m_heData[ m_heData[it.m_cur->next].twin];
		Eigen::Vector3d v2 = m_vertices[he.vert] - m_vertices[vertex];

		sum += this->angle_two_vector(v1, v2);

	} while(advance_iterator(it));

	//return sum / (2 * PI);
	if (sum < 2 * PI)
		return (1 - sum / (2 * PI));
	else
		return sum / (2 * PI);


	//return fabs(sum - 2 * PI);

	// for a better error
	double _sum(0);
	vvert_iterator iter;
	if (!this->init_iterator(iter, vertex))
		std::cerr << "error because of float point" << std::endl;

	do {

		std::size_t v1 = iter.m_cur->vert;
		std::size_t v2 = m_heData[iter.m_cur->twin].vert;
		_sum += get_quadric_error( *this, v1, v2 );
	} while(advance_iterator(iter));

	//if (m_numVertRemoved < 2 * m_vertData.size() / 3)
	//	return fabs(sum - 2 * PI);
	//else
	return 	_sum;
}

double EditMesh::get_angle_error(std::size_t vertex)
{
	double sum(0);
	Eigen::Vector3d vBarycent(0, 0, 0);

	// Boundary handling
	if (this->isBoundaryVert(vertex))
	{
		//std::cerr << "error because of boundary point" << std::endl;
		return std::numeric_limits<double>::max();
	}

	// initial the iterator visit the one ring of "vertex"
	vvert_iterator it;
	if (!this->init_iterator(it, vertex))
		std::cerr << "error because of float point" << std::endl;

	// sum the 1-ring angles	
	do {
		//double weight(this->get_mean_value_weight(it));
		double weight = 1;

		vBarycent += weight * m_vertices[it.m_cur->vert];

		sum += weight;

	} while(advance_iterator(it));

	return (m_vertices[vertex]  - vBarycent / sum).norm();	

	// for a better error
	double _sum(0);
	vvert_iterator iter;
	if (!this->init_iterator(iter, vertex))
		std::cerr << "error because of float point" << std::endl;

	do {

		std::size_t v1 = iter.m_cur->vert;
		std::size_t v2 = m_heData[iter.m_cur->twin].vert;
		_sum += get_quadric_error( *this, v1, v2 );
	} while(advance_iterator(iter));

	//if (m_numVertRemoved < 2 * m_vertData.size() / 3)
	//	return (m_vertices[vertex]  - vBarycent / sum).norm();	
	//else
	return 	_sum;
}

//determine which err to use
void EditMesh::Errormetric_simplify(bool mode)
{
	m_simplifyMetric = mode;

	// Initialization
	if (m_simplifyError.size() != m_vertData.size())
	{
		m_simplifyError.resize(m_vertData.size());
		for (std::size_t i(0); i < m_vertData.size(); ++i)
		{
			//gaussian error - point to plane
			if (m_simplifyMetric)
				m_simplifyError[i] = this->get_Gauss_error(i);
			//sum of angles error
			else
				m_simplifyError[i] = this->get_angle_error(i);
		}
	}

	// find first mins for coloring
	// first mins for coloring
	std::multimap<double, size_t> firstMins;	
	m_numColored = (m_numColored > m_vertData.size() - m_numVertRemoved) ?  m_vertData.size() - m_numVertRemoved : m_numColored;
	for (size_t i(0); i < m_numColored; ++i)
		firstMins.insert(std::make_pair(numeric_limits<double>::max(), HOLE_INDEX));

	for (size_t i(0), iEnd(m_simplifyError.size()); i < iEnd; ++i)
	{
		if (m_simplifyError[i] < firstMins.rbegin()->first)
		{
			auto iter = firstMins.end();
			firstMins.erase(--iter);
			firstMins.insert(make_pair(m_simplifyError[i], i));
		}
	}

	vertVecColor.clear();
	for (auto iter(firstMins.begin()); iter != firstMins.end(); ++iter)
	{
		if (iter->first != numeric_limits<double>::max())
			vertVecColor.push_back(iter->second);
	}
}

double EditMesh::angle_two_vector(Eigen::Vector3d v1, Eigen::Vector3d v2)
{
	v1.normalize();
	v2.normalize();
	double angle_cos = v1.dot(v2);

	if (fabs(angle_cos + 1) < 0.01)
		return PI;
	if (fabs(angle_cos - 1) < 0.01)
		return 0;

	return acos(angle_cos);
}

/*========================================
* simplification scheme algorithms 
*========================================*/
//set number of verts to be removed
void EditMesh::num_vert_remove(std::size_t num) // num: the number of vertice to be removed
{
	m_numVert = num;
}

void EditMesh::oneStep_simplify()
{
	//m_numVertRemoved : number of removed verts

	m_numVert = (m_numVert > m_vertData.size() - m_numVertRemoved) ?  m_vertData.size() - m_numVertRemoved : m_numVert;
	size_t _m_numVertRemoved(m_numVertRemoved);

	for (std::size_t i(0); i < m_numVert; ++i)
	{
		if (!this->remove_one_vertex())
			break;
	}
		
	if (_m_numVertRemoved != m_numVertRemoved)
	{
		std::cout <<m_numVertRemoved<<" removed vertices."<<endl;
		std::cout <<m_vertData.size() - m_numVertRemoved<<" remaining vertices."<<endl;
		std::cout <<"-----------------------------------"<<endl;
	}
}
bool EditMesh::remove_one_vertex()
{
	// Initialization
	if (m_simplifyError.size() != m_vertData.size())
	{
		m_simplifyError.resize(m_vertData.size());
		//init here
		//decide which error metric to use
		m_simplifyMetric = true;
		for (std::size_t i(0); i < m_vertData.size(); ++i)
		{
			/*if(m_vertData[i] == HOLE_INDEX)
			{
				m_simplifyError[i] = std::numeric_limits<double>::max();
				continue;
			}
*/
			if (m_simplifyMetric)
				m_simplifyError[i] = this->get_Gauss_error(i);
			else
				m_simplifyError[i] = this->get_angle_error(i);
		}
	}
	
	// find the minimal error vertex and remove it
	double minError(0);
	bool isEligibleRomove;
	std::vector<std::size_t> neighborVert;
	do{
		// Find the lowest error vertex.
		std::size_t vertIndex = HOLE_INDEX;
		auto iterIndex(std::min_element(m_simplifyError.begin(), m_simplifyError.end()));
		if (*iterIndex == std::numeric_limits<double>::max())
		{
			// can not be simplified anymore
			std::cerr << "ERROR! Can't simplify anymore!" << std::endl;
			vertVecColor.clear();
			return false;
			
		}
		else
		{
			vertIndex = std::distance(m_simplifyError.begin(), iterIndex);
			minError = *iterIndex;
		}

		// store the 1-ring vertex while our vert indices are still valid
		//i.e before deleting anything
		vvert_iterator it;		// initial the iterator visiting the one ring of "vertex"
		neighborVert.clear();
		if (!this->init_iterator(it,vertIndex))
		{
			std::cerr << "error because of float point" << std::endl;
			return false;
		}
		do {
			neighborVert.push_back(it.m_cur->vert);			
		} while(advance_iterator(it));
		
		// preform removing
		isEligibleRomove = this->vetex_removal_algorithm(vertIndex);
		//set the error of the removed vert to max
		m_simplifyError[vertIndex] = std::numeric_limits<double>::max();

		// store the deleted information for backward adding
		if (isEligibleRomove)
		{
			// store the removed vertex (at the end), and its 1-ring in a clockwise order
			neighborVert.push_back(vertIndex);	
			removed_verts.push(neighborVert);
			
			neighborVert.pop_back();
		}

	}while(!isEligibleRomove);

	++m_numVertRemoved;

	// Update the error values for all the old 1-ring vertices.
	minError /= neighborVert.size();
	for (auto iter(neighborVert.begin()); iter != neighborVert.end(); ++iter)
	{
		if (m_simplifyMetric)
			m_simplifyError[*iter] = minError + this->get_Gauss_error(*iter);
			
		else
			m_simplifyError[*iter] = minError + this->get_angle_error(*iter);
	}

	// first mins for coloring
	std::multimap<double, size_t> firstMins;	
	m_numColored = (m_numColored > m_vertData.size() - m_numVertRemoved) ?  m_vertData.size() - m_numVertRemoved : m_numColored;
	for (size_t i(0); i < m_numColored; ++i)
		firstMins.insert(std::make_pair(numeric_limits<double>::max(), HOLE_INDEX));
	
	for (size_t i(0), iEnd(m_simplifyError.size()); i < iEnd; ++i)
	{
		if (m_simplifyError[i] < firstMins.rbegin()->first)
		{
			auto iter = firstMins.end();
			firstMins.erase(--iter);
			firstMins.insert(make_pair(m_simplifyError[i], i));
		}
	}

	vertVecColor.clear();
	for (auto iter(firstMins.begin()); iter != firstMins.end(); ++iter)
	{
		if (iter->first != numeric_limits<double>::max())
		{
			vertVecColor.push_back(iter->second);
		}

	}

	this->edit_count_updater();

	return true;
}
/*========================================
* Checking if the vert can be removed;
*========================================*/
/*
1-check if vert is a boundry vert
2-check the change in normals
3-checks done in collapse he
4-check to ensure the shape remains manifold after removing the vert
*/
bool EditMesh::vetex_removal_algorithm(std::size_t vertIndex)
{
	//1 -------- Check if vertIndex is a boundary vertex
	//check by checking if the he is a boundry he
	if(isBoundaryVert(vertIndex))
		return false;
	/*========================================
	* Retriangulate the hole (fan)
	*========================================*/
	//find the start vertex which has the largest angle
	// find the best vertex from which fan trangulation strats 
	size_t bestStartHE;
	double minCost(std::numeric_limits<double>::max());
	vvert_iterator it; // initial the iterator visiting the one ring of "vertex"
	if (!this->init_iterator(it, vertIndex))
		std::cerr << "error because of float point" << std::endl;

	do {

		std::size_t v1 = it.m_cur->vert;
		std::size_t v2 = m_heData[it.m_cur->twin].vert;

		if (minCost > get_quadric_error( *this, v1, v2 ))
		bestStartHE = m_heData[it.m_cur->twin].twin;
	} while(advance_iterator(it));

	// check the topology eligibility
	//check the change in the normals
	it.m_cur = it.m_end = &m_heData[bestStartHE];
	while (advance_iterator(it))
	{
		// degree calculation
		//++degree;
		
		Eigen::Vector3d v1( m_vertices[m_heData[m_heData[it.m_cur->next].next].vert] - m_vertices[it.m_end->vert]);
		Eigen::Vector3d v2(m_vertices[it.m_cur->vert] - m_vertices[it.m_end->vert]);
		v1.normalize();
		v2.normalize();
		//double a = (v1.cross(v2)).dot(this->get_vnormal(vertIndex));
		if ((v1.cross(v2)).dot(this->get_vnormal(vertIndex)) < -0.5) // non convex condition
			return false;
	}	

	//2	----------From hedge collapse
	//size_t he(m_heData[m_vertData[vertIndex]].twin);
	size_t he(bestStartHE);
	//Best half-edge to start retriangulating from
	const half_edge& heBase = m_heData[he];
	const half_edge& heTwin = m_heData[heBase.twin];

	// We are going to delete the faces on either side of the chosen edge, so we need to delete 3 half_edges and patch up the twin links on the 4
	// bordering edges.
	std::size_t heBorder[4];
	heBorder[0] = m_heData[ heBase.next ].twin;
	heBorder[1] = m_heData[ m_heData[ heBase.next ].next ].twin;
	heBorder[2] = m_heData[ m_heData[ heTwin.next ].next ].twin;
	heBorder[3] = m_heData[ heTwin.next ].twin;

	// TODO: Relax this assertion. We should be able to collapse a spike jutting into a hole.
	assert( ( m_heData[ heBorder[0] ].face != HOLE_INDEX || m_heData[ heBorder[1] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );
	assert( ( m_heData[ heBorder[2] ].face != HOLE_INDEX || m_heData[ heBorder[3] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );

	// Check if we can actually collapse. This checks for a degree 3 vertex not on the edge we are collapsing.
	if( m_heData[ m_heData[ m_heData[ heBorder[1] ].next ].twin ].next == heBorder[0] ||
		m_heData[ m_heData[ m_heData[ heBorder[2] ].next ].twin ].next == heBorder[3] )
		return false;

	//3------- check non-manifold
	size_t degree(3);
	it.m_cur = &m_heData[bestStartHE];
	it.m_end = &m_heData[m_heData[m_heData[m_heData[bestStartHE].twin].next].next];
	advance_iterator(it);
	while (advance_iterator(it))
	{
		// the degree of vertex
		++degree;
		size_t v1(it.m_cur->vert);

		vvert_iterator itTemp; // initiate the iterator visiting the one ring of "vertex"
		if (!this->init_iterator(itTemp, m_heData[bestStartHE].vert))
			std::cerr << "error " << std::endl;
		do {
				size_t v2(itTemp.m_cur->vert);
				if (v1 == v2)	 // v1 and v2 should not be the same
					return false;
		} while(advance_iterator(itTemp));
	}	
	//if (degree == 3) return false; // checks for a degree 3 vertex not on the edge we are collapsing.


	// Capture the indices of things (2 faces & 6 half-edges) we want to delete.
	std::size_t fToDelete[] = { heBase.face, heTwin.face };
	std::size_t heToDelete[] = { he, heBase.next, m_heData[ heBase.next ].next, heBase.twin, heTwin.next, m_heData[ heTwin.next ].next };

	// We may also need to fix the vertex->half_edge link for the verts using these faces. 
	//There are technically 4, but we only update the 3 that are not going to be deleted.
	std::size_t verts[] = { this->prev( heBase ).vert, heBase.vert, this->prev( heTwin ).vert };

	// Adjust all the twin's 1-ring to link to the vertex we are not going to delete.
	std::size_t heIt = this->twin(this->next(heBase)).next;
	std::size_t heEnd = heBase.twin;
	//CW itteration of Twin hes and ensuring their verts is the middle vert
	for( ; heIt != heEnd; heIt = this->twin( m_heData[heIt] ).next ){
		assert( m_heData[heIt].vert == heTwin.vert );

		// now set their verts to be the heBase vert, so we can delete this one.
		m_heData[heIt].vert = heBase.vert;
	}

	// Fix the vert associations if required, picking a non-hole face.
	if( m_vertData[ verts[0] ] == m_heData[ heBorder[1] ].twin )
		m_vertData[ verts[0] ] = (m_heData[ heBorder[0] ].face != HOLE_INDEX) ? heBorder[0] : m_heData[ heBorder[1] ].next;
	if( m_vertData[ verts[1] ] == he || m_vertData[ verts[1] ] == heTwin.next )
		m_vertData[ verts[1] ] = (m_heData[ heBorder[1] ].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
	if( m_vertData[ verts[2] ] == m_heData[ heBorder[2] ].twin )
		m_vertData[ verts[2] ] = (m_heData[ heBorder[3] ].face != HOLE_INDEX) ? heBorder[3] : m_heData[ heBorder[2] ].next;

	// Collapse the two triangles bordering our chosen half-edge by connecting the opposite edges together.
	m_heData[ heBorder[0] ].twin = heBorder[1];
	m_heData[ heBorder[1] ].twin = heBorder[0];
	m_heData[ heBorder[2] ].twin = heBorder[3];
	m_heData[ heBorder[3] ].twin = heBorder[2];

	// "Delete" the vertex
	m_vertData[vertIndex] = HOLE_INDEX;

	// Have to delete the faces in the proper order.
	if( fToDelete[0] < fToDelete[1] )
		std::swap( fToDelete[0], fToDelete[1] );

	this->delete_half_edges_impl( heToDelete );
	detail::delete_face( m_faceData, m_heData, fToDelete[0] );
	detail::delete_face( m_faceData, m_heData, fToDelete[1] );

	for (size_t i(0), iEnd(m_heData.size()); i < iEnd; ++i	)
	{
		half_edge hfBase = m_heData[i];
		half_edge hfNext = m_heData[hfBase.next];
		half_edge hfNextNext = m_heData[hfNext.next];

		/*if( hfNextNext.next != i ){
			cout<<"Half edges don't form a loop"<<endl;

		}*/

		//assert( hfNextNext.next == i && "Half edges don't form a loop" );
	
		assert( m_heData[hfBase.twin].twin == i && "Twin half edges are not mutual" );
		assert( m_heData[hfNext.twin].twin == hfBase.next && "Twin half edges are not mutual");
		assert( m_heData[hfNextNext.twin].twin == hfNext.next && "Twin half edges are not mutual");

		assert( hfBase.face == hfNext.face && hfBase.face == hfNextNext.face && "Half edges are not in the same face" );
		assert( m_faceData[hfBase.face] ==i || m_faceData[hfBase.face] == hfBase.next || m_faceData[hfBase.face] == hfNext.next && "Half edges are not in the same face" );
	
		bool isVertPoint2HE(false);
		vvert_iterator itTemp;
		itTemp.m_cur = itTemp.m_end = &m_heData[hfBase.twin];
		do 
		{
			if (itTemp.m_cur->twin == m_vertData[hfBase.vert])
				isVertPoint2HE = true;

			assert (m_heData[itTemp.m_cur->twin].vert == hfBase.vert && "1-ring half edges don't share the same vertex");
		} while (advance_iterator(itTemp));

		assert (isVertPoint2HE && "Center vertex points to other bunch of half edges");	
	}

	return true;
}

/*========================================
* backwards simplification scheme algorithms
*(buggy)num_vert_remove
*========================================*/

void EditMesh::mesh_backward_simplify()
{
	// Check if the to be added number is greater than removed vertices
	m_numVert = (m_numVert > m_numVertRemoved) ? m_numVertRemoved : m_numVert;

	for (std::size_t i(0); i < m_numVert; ++i)
	{
		if (!this->add_back_one_vertex())
			break;
	}

	if (m_numVertRemoved == 0)
		std::cout << "Back to the original mesh - no more vertices to be added" << std::endl;
	else
		std::cout << m_numVert <<" vertices have been added and total number of vertices in the mesh are: " 
		<< m_vertData.size() - m_numVertRemoved << std::endl;
}

bool EditMesh::add_back_one_vertex()
{
	if (removed_verts.empty())
		
			return false;

	// Get the centric vertex and its 1-ring
	size_t centVert(removed_verts.top().back());
	removed_verts.top().pop_back();
	vector<size_t> neighborVert(removed_verts.top());
	removed_verts.pop();

	// find the bordering halfedges
	vector<half_edge> neighborHE;
	for (size_t i(0), iEnd(neighborVert.size()); i < iEnd; ++i)
	{
		size_t vFrom(neighborVert[i]);
		size_t vTo(neighborVert[(i + 1) % iEnd]);

		half_edge* hePtr(this->find_edge(vFrom, vTo));
		if (hePtr == NULL)
		{
			// can not find the dege, this is impossible.
			std::cerr << "There is error in the half edge data!" << std::endl;
			return false;
		}
		else
			neighborHE.push_back(*hePtr);
	}

	// half edges and faces to be deleted
	// Capture the indices of things (2 faces & 6 half-edges) we want to delete.
	set<size_t> fToDelete, heToDelete;
	for (auto iter(neighborHE.begin()); iter != neighborHE.end(); ++iter)
	{
		// half edges
		heToDelete.insert(iter->twin);
		heToDelete.insert(m_heData[iter->twin].next);
		heToDelete.insert(m_heData[m_heData[iter->twin].next].next);

		// faces
		fToDelete.insert(m_heData[iter->twin].face);
	}

	// update half edges, faces and vertex
	size_t heNum(m_heData.size());
	size_t faceNum(m_faceData.size());
	m_heData.resize(heNum + 3 * neighborHE.size());
	m_faceData.resize(faceNum + neighborHE.size());

	for (size_t i(0), iEnd(neighborHE.size()); i < iEnd; ++i)
	{
		// attach the from vertex to the half edge
		size_t indexCurHE(m_heData[neighborHE[i].twin].twin);
		m_vertData[neighborHE[i].vert] = indexCurHE;
		
		size_t twin(heNum + 3 * i + 0);
		size_t next(heNum + 3 * i + 1);
		size_t prev(heNum + 3 * i + 2);

		size_t twinFromNextFace(heNum + 3 * ((i + 1) % neighborHE.size()) + 1);
		size_t twinFromPrevFace(heNum + 3 * ((i + neighborHE.size() - 1) % neighborHE.size()) + 2);

		// new half edges
		detail::init( m_heData[twin], next, indexCurHE, neighborVert[(i + 1) % iEnd], faceNum + i);
		detail::init( m_heData[next], prev, twinFromPrevFace, neighborHE[i].vert, faceNum + i);
		detail::init( m_heData[prev], twin, twinFromNextFace, centVert, faceNum + i);

		m_heData[indexCurHE].twin = twin;

		// new face
		m_faceData[faceNum + i] = twin;

		// add the half edge for center vertex
		m_vertData[centVert] = prev;
	}

	for (auto iter(heToDelete.rbegin()); iter != heToDelete.rend(); ++iter)
	{
		double i = *iter;
		this->delete_half_edge_impl(*iter);
	}

	for (auto iter(fToDelete.rbegin()); iter != fToDelete.rend(); ++iter)
		detail::delete_face( m_faceData, m_heData, *iter);

	if (m_numVertRemoved == 0)
		return false;
	else
		--m_numVertRemoved;

	// Update the error values for all the old 1-ring vertices.
	for (auto iter(neighborVert.begin()); iter != neighborVert.end(); ++iter)
	{
		if (m_simplifyMetric)
			m_simplifyError[*iter] = this->get_Gauss_error(*iter);
		else
			m_simplifyError[*iter] = this->get_angle_error(*iter);
	}

	m_simplifyError[centVert] = this->get_angle_error(centVert);

	// first mins for coloring
	std::multimap<double, size_t> firstMins;	
	m_numColored = (m_numColored > m_vertData.size() - m_numVertRemoved) ?  m_vertData.size() - m_numVertRemoved : m_numColored;
	for (size_t i(0); i < m_numColored; ++i)
		firstMins.insert(std::make_pair(numeric_limits<double>::max(), HOLE_INDEX));

	for (size_t i(0), iEnd(m_simplifyError.size()); i < iEnd; ++i)
	{
		if (m_simplifyError[i] < firstMins.rbegin()->first)
		{
			auto iter = firstMins.end();
			firstMins.erase(--iter);
			firstMins.insert(make_pair(m_simplifyError[i], i));
		}
	}

	vertVecColor.clear();
	for (auto iter(firstMins.begin()); iter != firstMins.end(); ++iter)
	{
		if (iter->first != numeric_limits<double>::max())
		{
			vertVecColor.push_back(iter->second);
		}

	}

	this->edit_count_updater();

	return true;
}

void EditMesh::simplificationScheme(){

	size_t numEdges = m_heData.size()/2;
	size_t numVerts = m_vertData.size();
	std::vector<bool> heChecker(numEdges,false) ;

	
		
	std::vector< size_t > IHeVec;
	std::vector< double > edge_lenVec;
	std::vector< Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > edgeVec;
	std::vector< size_t > NumHeInOneRing;
	// angle between to edge lengths 
	std::vector< double > angleVec; 
	double angleSum = 0;
	

	// allocate memory
	vertErr* localErrArr = static_cast<vertErr*>( ::operator new ( sizeof vertErr * UnsortedErr.size()) );
	//std::vector<vertErr>localErrUnsortedArr;
	//vertErr* localErrUnsortedArr = static_cast<vertErr*>( ::operator new ( sizeof vertErr * m_vertData.size()) );

	for (size_t i=0; i<UnsortedErr.size(); i++)
	{
		
		//if (localErrUnsortedArr[i]._vertI!=HOLE_INDEX){
		   /*========================================
			* Simple Error
			*========================================*/
			//Store half edges originating form  the specified vertex.
			cout<<UnsortedErr[i]._vertI<<endl;
			IHeVec = IHeInRing(UnsortedErr[i]._vertI);
			
			/*========================================
			* sum of incident angles.
			*========================================*/
			for (size_t j=0; j<IHeVec.size(); j++)
			{
				const half_edge *itCur = &m_heData[IHeVec[j]];
				const half_edge *itTwin = &m_heData[ itCur->twin ];
				Eigen::Vector3d a = this->get_vertex( itTwin->vert );
				//original vertex
				Eigen::Vector3d b = this->get_vertex( itCur->vert );

				assert( ( itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX ) && "Invalid mesh: edge with no face on either side" );

				if( itCur->face != HOLE_INDEX ){
					Eigen::Vector3d c = this->get_vertex( this->prev( *itCur ).vert );
	
					Eigen::Vector3d e0 = a - b;
					Eigen::Vector3d e1 = c - b;
					edgeVec.push_back(e0);

					// edge lengths of the two triangles, used to determine the 2d representation
					double    e0_len = e0.norm();
					double    e1_len = e1.norm();
					edge_lenVec.push_back(e0_len);

					// We use the dot product to get cos 
					double  angle = static_cast<double>(std::acos(( e0.dot(e1) ) / (e0_len*e1_len)));
					//cout << "angle"<< j+1 <<"=\n" << angle << std::endl;
					angleVec.push_back(angle);
				
					angleSum += angle;
				
			}
			
				}
			/*double localErr = abs(angleSum/(2*M_PI)-1);*/

			/*========================================
			* The more complicated error fnc
			*========================================*/
			//get the border verts
			vector<size_t>bVertsForErr;
			for (size_t j=0; j<IHeVec.size(); j++)
			{
				const half_edge *itCur = &m_heData[IHeVec[j]];
				const half_edge *itTwin = &m_heData[ itCur->twin ];
				bVertsForErr.push_back(itTwin->vert);
			}
			//average the border verts to get the base point
			Eigen::Vector3d bVertsSum;
			Eigen::Vector3d planeBasePoint;
			for(size_t j=0;j<bVertsForErr.size();j++){
				bVertsSum += m_vertices[bVertsForErr[j]];
			}
			planeBasePoint = bVertsSum/bVertsForErr.size();

			Eigen::Vector3d posToPlane = m_vertices[UnsortedErr[i]._vertI] - planeBasePoint;
			double distToPlane = posToPlane.norm();
			/*========================================*/

			double localErr = distToPlane;
			cout<<"local Err "<<localErr<<endl;
			localErrArr[i] = vertErr( localErr,UnsortedErr[i]._vertI  );
			//cout << "sum of the agles =\n" << angleSum << std::endl;

		
			NumHeInOneRing.push_back(IHeVec.size());
			angleSum = 0;
			IHeVec.clear();
			angleVec.clear();
			cout<<"done! round "<<i<<endl;
	}
	
	 // Output the heap after each element that we add
	for (size_t i = 0; i < UnsortedErr.size(); ++i) {
		UnsortedErr[i]._vertI = localErrArr[i]._vertI;
		UnsortedErr[i]._vertError = UnsortedErr[i]._vertError + 
									localErrArr[i]._vertError;

		AddElement(localErrArr[i], localErrArr, i);
	}
	
	// Output the heap after each element that we remove
	for (int i = 0; i < UnsortedErr.size(); ++i) {
		//first element has the smallest error
		RemoveRoot(localErrArr, UnsortedErr.size() - i);
		
	}

	
	/*========================================
	* Checking if the vert can be removed;
	* if so, remove the vert and other 
	* necessary faces and hes
	*========================================*/
	//checking 
	//if vToDelete ends up being a HOLE_INDEX, there is no more vertices to remove
	
	size_t* vToDelete = NULL;
	//vertices on the border
	vector <size_t> borderverts; 
	vector<size_t>heBorder;
	double localErr = 0;
	for (size_t i=0; i<UnsortedErr.size(); i++){
		borderverts.clear();
		//check the first condition
		localErr = localErrArr[i]._vertError;
		vToDelete = &localErrArr[i]._vertI;
		// get the bordering half-edges.
		heBorder = borderHes(*vToDelete);
		cout<<heBorder.size()<<endl;
		
		if(localErrArr[i]._vertI==HOLE_INDEX)
			continue;


	    //Get the border vertices
		vector <size_t> orgHeVec = IHeInRing(localErrArr[i]._vertI);
		
		//cout<<"The border verts are "<<endl;
		for (size_t j=0; j<orgHeVec.size(); j++)
		{
			const half_edge *itCur = &m_heData[orgHeVec[j]];
			const half_edge *itTwin = &m_heData[ itCur->twin ];
			borderverts.push_back(itTwin->vert);
			
			/*cout<<borderverts[j];
			cout<<endl;*/
		}
		
		//select one of the border vertices to retriangulate
		//???Do we need to do any checks on this vert?
		size_t vBorder = borderverts[0];


		//size_t vertBorder = m_heData[m_heData[m_vertData[*vToDelete]].twin].vert;
		half_edge* he = find_edge(vBorder,*vToDelete);
		const half_edge& heBase = * he;
		const half_edge& heTwin = m_heData[heBase.twin];

		// We are going to delete the faces on either side of the chosen edge, so we need to delete 3 half_edges and patch up the twin links on the 4
		// bordering edges.
		std::size_t heB[4];
		heB[0] = m_heData[ heBase.next ].twin;
		heB[1] = m_heData[ m_heData[ heBase.next ].next ].twin;
		heB[2] = m_heData[ m_heData[ heTwin.next ].next ].twin;
		heB[3] = m_heData[ heTwin.next ].twin;

		// Check if we can actually collapse. This checks for a degree 3 vertex at the vertices not on the edge we are collapsing.
		if( m_heData[ m_heData[ m_heData[ heB[1] ].next ].twin ].next == heB[0] ){
			*vToDelete = HOLE_INDEX;
			continue;
		}
		if( m_heData[ m_heData[ m_heData[ heB[2] ].next ].twin ].next == heB[3] ){
			*vToDelete = HOLE_INDEX;
			continue;
		}



		//// Check if we can actually collapse. This checks for a degree 3 vertex at the vertices not on the edge we are collapsing.
		//if( m_heData[ m_heData[ m_heData[ heBorder[1] ].next ].twin ].next == heBorder[0] ){
		//	*vToDelete = HOLE_INDEX;
		//	continue;
		//}
		//if( m_heData[ m_heData[ m_heData[ heBorder[2] ].next ].twin ].next == heBorder[3] ){
		//	*vToDelete = HOLE_INDEX;
		//	continue;
		//}

		
		//check if the mesh becomes non-manifold
		////Get the border vertices
		//vector <size_t> orgHeVec = IHeInRing(localErrArr[i]._vertI);
		//
		////cout<<"The border verts are "<<endl;
		//for (size_t j=0; j<orgHeVec.size(); j++)
		//{
		//	const half_edge *itCur = &m_heData[orgHeVec[j]];
		//	const half_edge *itTwin = &m_heData[ itCur->twin ];
		//	borderverts.push_back(itTwin->vert);
		//	
		//	/*cout<<borderverts[j];
		//	cout<<endl;*/
		//}
		//
		////select one of the border vertices to retriangulate
		////???Do we need to do any checks on this vert?
		//size_t vBorder = borderverts[0];
		/*cout<<"The selected border vert for retriangulation is "<<vBorder;
		cout<<endl;*/

		//borderverts[1] = HOLE_INDEX;
		//borderverts[orgHeVec.size()-1] = HOLE_INDEX;
		/*cout<<"The excluded border verts are"<<borderverts[1]
			<<" and "<<borderverts[orgHeVec.size()-1];
		cout<<endl;*/

		vector <size_t> orgHeVecFromVBorder = IHeInRing(vBorder);
		vector <size_t> NVertsofBorderVert; 
		/*cout<<"size of the orgHeVecFromVBorder is "<<orgHeVecFromVBorder.size();
		cout<<endl;*/

		/*cout<<"neighbor verts of border vert are";
		cout<<endl;*/
		for (size_t j=0; j<orgHeVecFromVBorder.size(); j++)
		{
			const half_edge *itCur = &m_heData[orgHeVecFromVBorder[j]];
			const half_edge *itTwin = &m_heData[ itCur->twin ];
			NVertsofBorderVert.push_back(itTwin->vert);
			/*cout<<NVertsofBorderVert[j];
			cout<<endl;*/
		}
		//loop breaks if the vert to remove becomes a whole index
		for (size_t j=2; j<borderverts.size()-1; j++){
			if (*vToDelete == HOLE_INDEX)
				break;
			for (int k=0; k< NVertsofBorderVert.size(); k++){
				if (borderverts[j] == NVertsofBorderVert[k]){
					cout<<"borderverts "<<borderverts[j]<<endl;
					cout<<"NVertsofBorderVert "<<borderverts[k]<<endl;
					*vToDelete = HOLE_INDEX;
					break;
				}
			}
		}
		if (*vToDelete == HOLE_INDEX)
			continue;
		
		break;
	}
	//needs to change later
	if (*vToDelete==HOLE_INDEX){
			cout<<"no available vertices to remove"<<endl;
			return;
	}
	else
			cout<<"vert to remove is : "<<*vToDelete<<endl;
	
	if (*vToDelete!=HOLE_INDEX){

		cout<<"borderverts.size() "<<borderverts.size()<<endl;
		cout<<"UnsortedErr.size()"<<UnsortedErr.size()<<endl;

		//update the error on the neighboring verts of the romoved
		//vert for the next round of itteration
		double updatedErr = localErr/borderverts.size();
		for(size_t i=0;i<borderverts.size();i++){
			for(size_t j=0;j<UnsortedErr.size();j++){
				if (UnsortedErr[j]._vertI == borderverts[i]){
					UnsortedErr[j]._vertError = UnsortedErr[j]._vertError + updatedErr;
					break;
				}

			}
			/*UnsortedErr[borderverts[i]]._vertError =
				UnsortedErr[borderverts[i]]._vertError + updatedErr;*/
				//cout<<UnsortedErr[borderverts[i]]._vertError<<endl;
		}
		cout<<"Line 2095 Done!"<<endl;
		//save necessary stuff before deleting anything
		vector <half_edge> inHeBorderObj = InnerBorderHes( *vToDelete);

		//cout<<"inHeBorderObj size is "<<inHeBorderObj.size()<<endl;
		//cout<<"m_faceData.size() before deleting "<<m_faceData.size()<<endl;
		
		//---------deleting all the faces---------------
	
		vector<size_t>hesToDelete = IHeInRing(*vToDelete);
		//vector of neighbouring faces around the vert 
		vector<size_t> fNeighbor ;
		
		//cout<<m_faceData.size()<<endl;
		for (size_t i=0;i<hesToDelete.size();i++){
			fNeighbor.push_back(m_heData[hesToDelete[i]].face);
			//cout<<fNeighbor[i]<<endl;
			//detail::delete_face( m_faceData, m_heData, fNeighbor[i]);
			//cout<<m_heData[hesToDelete[i]].face<<endl;
		}
		deleteUmbrellaFaces(*vToDelete);

		//----------deleting hedges----------------
		/*vector<size_t>hesToDelete = IHeInRing(*vToDelete);*/

		//NOTE: for it to work, half-edges must be removed in a descending order
		/*for(size_t i=hesToDelete.size()-1;i>-1;i--)
			this->delete_half_edge_impl(hesToDelete[i]);*/

		//-------deleting the vertix---------
		/*for(size_t i=0;i<UnsortedErr.size();i++){
			cout<<UnsortedErr[i]._vertI<<" "<<UnsortedErr[i]._vertError<<endl;
	}*/
		for(size_t j=0;j<UnsortedErr.size();j++){
			if (UnsortedErr[j]._vertI == *vToDelete){
				UnsortedErr[j]._vertI =  HOLE_INDEX;
				UnsortedErr.erase(UnsortedErr.begin()+j);
				break;
			}
		}
		/*UnsortedErr[*vToDelete]._vertI =  HOLE_INDEX;
		UnsortedErr.erase(UnsortedErr.begin()+*vToDelete);*/
		//UnsortedErr[*vToDelete]._vertError = 100000
		/*for(size_t i=0;i<UnsortedErr.size();i++)
			cout<<UnsortedErr[i]._vertI<<" "<<UnsortedErr[i]._vertError<<endl;*/
	/*========================================
	* Retriangulate the hole
	*========================================*/
	//size_t newFaceI = this->add_face(firstHe->vert, secondHe->vert, thirdHe.vert);

	//we need to create the new hes
	vector <half_edge> newHes;
	vector <size_t> InewHes;
	vector <half_edge> newTwinHes;
	vector <size_t> InewTwinHes;
	vector <size_t>bordervertsmod;
	vector <half_edge*>HeBorder;
	vector <size_t>inHeBorder;
	
	vector <size_t> newFacesI;

	//recreating the faces
	for(size_t i=0;i<fNeighbor.size()-2;i++){
		newFacesI.push_back(m_faceData.size()+i);
		//cout<<newFacesI[i]<<endl;
	}
	cout<<"line 3182 Done"<<endl;

	//finding new indecis of border hes after deleting

	for (size_t i=0;i<borderverts.size()-1;i++){
		HeBorder.push_back(find_edge( borderverts[i], borderverts[i+1] ));		
		//cout<<HeBorder[i]->face<<endl;
	}
	HeBorder.push_back(find_edge( borderverts[borderverts.size()-1], borderverts[0] ));
	
	for (size_t i=0;i<HeBorder.size();i++){
		inHeBorder.push_back(m_heData[HeBorder[i]->twin].twin);
		//cout<<"inHeBorder are "<<inHeBorder[i]<<endl;
		//cout<<inHeBorder[i]<<endl;
	}

	//first face -> face[0]
	//cout<<m_heData[inHeBorder[0]].face<<endl;
	m_faceData.push_back(inHeBorder[0]);
	m_heData[inHeBorder[0]].face = newFacesI[0];
	m_heData[inHeBorder[0]].next = inHeBorder[1];
	

	//cout<<"B0 vert is: "<<m_heData[inHeBorder[0]].vert<<endl;
	
	for (size_t i =1; i<inHeBorder.size()-2; i++){
		//creating the new half edge and its twin
		newHes.push_back(half_edge());
		InewHes.push_back(m_heData.size());
		newTwinHes.push_back(half_edge());
		InewTwinHes.push_back(m_heData.size()+1);
		
		//inner border half-edge
		m_heData[inHeBorder[i]].face = newFacesI[i-1];
		m_heData[inHeBorder[i]].next = InewHes[i-1];

		//new half-edge
		newHes[i-1].vert = borderverts[i+1];
		/*cout<<"new half-edge vert: "<<endl;
		cout<<newHes[i-1].vert <<endl;
		cout<<m_heData[m_heData[inHeBorder[i]].twin].vert <<endl;*/

		newHes[i-1].face = newFacesI[i-1];
		newHes[i-1].twin = InewTwinHes[i-1];
		if(i==1)
			newHes[i-1].next = inHeBorder[i-1];
		else
			newHes[i-1].next = newHes[i-2].twin;

		//new twin half-edge
		m_faceData.push_back(InewTwinHes[i-1]);
		newTwinHes[i-1].vert = borderverts[0];
		newTwinHes[i-1].face = newFacesI[i];
		newTwinHes[i-1].twin = InewHes[i-1];
		newTwinHes[i-1].next = inHeBorder[i+1]; 
		
		m_heData.push_back(newHes[i-1]);
		m_heData.push_back(newTwinHes[i-1]);



	}
	//the last two inner border half-edges

	m_heData[inHeBorder[inHeBorder.size()-1]].face = newFacesI[newFacesI.size()-1];
	m_heData[inHeBorder[inHeBorder.size()-1]].next = InewTwinHes[InewTwinHes.size()-1];
	//cout<<"vert of  last bhe is: "<<m_heData[inHeBorder[inHeBorder.size()-1]].vert<<endl;

	m_heData[inHeBorder[inHeBorder.size()-2]].face = newFacesI[newFacesI.size()-1];
	m_heData[inHeBorder[inHeBorder.size()-2]].next = inHeBorder[inHeBorder.size()-1];
	//cout<<"vert of one to last bhe is: "<<m_heData[inHeBorder[inHeBorder.size()-2]].vert<<endl;
	
	//cout<<"m_faceData.size() after deleting"<<m_faceData.size()<<endl;

	
	/*cout<<"new hes twins verts are"<<endl;
	for (size_t i=0;i<newHes.size();i++)
		cout<<newTwinHes[i].vert<<endl;*/

	fNeighbor.clear();
	

	}
	cout<<"end"<<endl;
	delete[] localErrArr;
	++edit_count;
}


bool EditMesh::remove_vetex(std::size_t vertIndex)
{
	// Check if vertIndex is a boundary vertex
	if(isBoundaryVert(vertIndex))
		return false;

	//// fan shaped triangulation of hole: find the start vertex which has the largest angle
	// find the best vertex from which fan trangulation strats 
	size_t bestStartHE;
	double minCost(std::numeric_limits<double>::max());
	vvert_iterator it; // initial the iterator visiting the one ring of "vertex"
	if (!this->init_iterator(it, vertIndex))
		std::cerr << "error because of float point" << std::endl;

	do {

		std::size_t v1 = it.m_cur->vert;
		std::size_t v2 = m_heData[it.m_cur->twin].vert;

		if (minCost > get_quadric_error( *this, v1, v2 ))
		bestStartHE = m_heData[it.m_cur->twin].twin;
	} while(advance_iterator(it));

	// check the topology eligibility
	it.m_cur = it.m_end = &m_heData[bestStartHE];
	while (advance_iterator(it))
	{
		// degree calculation
		//++degree;
		
		Eigen::Vector3d v1( m_vertices[m_heData[m_heData[it.m_cur->next].next].vert] - m_vertices[it.m_end->vert]);
		Eigen::Vector3d v2(m_vertices[it.m_cur->vert] - m_vertices[it.m_end->vert]);
		v1.normalize();
		v2.normalize();
		//double a = (v1.cross(v2)).dot(this->get_vnormal(vertIndex));
		if ((v1.cross(v2)).dot(this->get_vnormal(vertIndex)) < -0.5) // non convex condition
			return false;
	}	

	//size_t he(m_heData[m_vertData[vertIndex]].twin);
	size_t he(bestStartHE);

	const half_edge& heBase = m_heData[he];
	const half_edge& heTwin = m_heData[heBase.twin];

	// We are going to delete the faces on either side of the chosen edge, so we need to delete 3 half_edges and patch up the twin links on the 4
	// bordering edges.
	std::size_t heBorder[4];
	heBorder[0] = m_heData[ heBase.next ].twin;
	heBorder[1] = m_heData[ m_heData[ heBase.next ].next ].twin;
	heBorder[2] = m_heData[ m_heData[ heTwin.next ].next ].twin;
	heBorder[3] = m_heData[ heTwin.next ].twin;

	// TODO: Relax this assertion. We should be able to collapse a spike jutting into a hole.
	assert( ( m_heData[ heBorder[0] ].face != HOLE_INDEX || m_heData[ heBorder[1] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );
	assert( ( m_heData[ heBorder[2] ].face != HOLE_INDEX || m_heData[ heBorder[3] ].face != HOLE_INDEX ) && "Cannot collapse an edge on a face with holes on either side." );

	// Check if we can actually collapse. This checks for a degree 3 vertex not on the edge we are collapsing.
	if( m_heData[ m_heData[ m_heData[ heBorder[1] ].next ].twin ].next == heBorder[0] ||
		m_heData[ m_heData[ m_heData[ heBorder[2] ].next ].twin ].next == heBorder[3] )
		return false;

	// check non-manifold
	size_t degree(3);
	it.m_cur = &m_heData[bestStartHE];
	it.m_end = &m_heData[m_heData[m_heData[m_heData[bestStartHE].twin].next].next];
	advance_iterator(it);
	while (advance_iterator(it))
	{
		// the degree of vertex
		++degree;
		size_t v1(it.m_cur->vert);

		vvert_iterator itTemp; // initial the iterator visiting the one ring of "vertex"
		if (!this->init_iterator(itTemp, m_heData[bestStartHE].vert))
			std::cerr << "error because of float point" << std::endl;
		do {
				size_t v2(itTemp.m_cur->vert);
				if (v1 == v2)	 // v1 and v2 should not be the same
					return false;
		} while(advance_iterator(itTemp));
	}	
	//if (degree == 3) return false; // checks for a degree 3 vertex not on the edge we are collapsing.


	// Capture the indices of things (2 faces & 6 half-edges) we want to delete.
	std::size_t fToDelete[] = { heBase.face, heTwin.face };
	std::size_t heToDelete[] = { he, heBase.next, m_heData[ heBase.next ].next, heBase.twin, heTwin.next, m_heData[ heTwin.next ].next };

	// We may also need to fix the vertex->half_edge link for the verts using these faces. There are technically 4, but we only update the 3 that are not going to be deleted.
	std::size_t verts[] = { this->prev( heBase ).vert, heBase.vert, this->prev( heTwin ).vert };

	// Adjust all the twin's 1-ring to link to the vertex we are not going to delete.
	std::size_t heIt = this->twin(this->next(heBase)).next;
	std::size_t heEnd = heBase.twin;
	for( ; heIt != heEnd; heIt = this->twin( m_heData[heIt] ).next ){
		assert( m_heData[heIt].vert == heTwin.vert );

		// Associate to the other vertex now, so we can delete this one.
		m_heData[heIt].vert = heBase.vert;
	}

	// Fix the vert associations if required, picking a non-hole face.
	if( m_vertData[ verts[0] ] == m_heData[ heBorder[1] ].twin )
		m_vertData[ verts[0] ] = (m_heData[ heBorder[0] ].face != HOLE_INDEX) ? heBorder[0] : m_heData[ heBorder[1] ].next;
	if( m_vertData[ verts[1] ] == he || m_vertData[ verts[1] ] == heTwin.next )
		m_vertData[ verts[1] ] = (m_heData[ heBorder[1] ].face != HOLE_INDEX) ? heBorder[1] : heBorder[2];
	if( m_vertData[ verts[2] ] == m_heData[ heBorder[2] ].twin )
		m_vertData[ verts[2] ] = (m_heData[ heBorder[3] ].face != HOLE_INDEX) ? heBorder[3] : m_heData[ heBorder[2] ].next;

	// Collapse the two triangles bordering our chosen half-edge by connecting the opposite edges together.
	m_heData[ heBorder[0] ].twin = heBorder[1];
	m_heData[ heBorder[1] ].twin = heBorder[0];
	m_heData[ heBorder[2] ].twin = heBorder[3];
	m_heData[ heBorder[3] ].twin = heBorder[2];

	// "Delete" the vertex
	m_vertData[vertIndex] = HOLE_INDEX;

	// Have to delete the faces in the proper order.
	if( fToDelete[0] < fToDelete[1] )
		std::swap( fToDelete[0], fToDelete[1] );

	this->delete_half_edges_impl( heToDelete );
	detail::delete_face( m_faceData, m_heData, fToDelete[0] );
	detail::delete_face( m_faceData, m_heData, fToDelete[1] );

	// check the data structure of mesh
	for (size_t i(0), iEnd(m_heData.size()); i < iEnd; ++i	)
	{
		half_edge hfBase = m_heData[i];
		half_edge hfNext = m_heData[hfBase.next];
		half_edge hfNextNext = m_heData[hfNext.next];

		assert( hfNextNext.next == i && "Half edges don't form a loop" );
		assert( m_heData[hfBase.twin].twin == i && "Twin half edges are not mutual" );
		assert( m_heData[hfNext.twin].twin == hfBase.next && "Twin half edges are not mutual");
		assert( m_heData[hfNextNext.twin].twin == hfNext.next && "Twin half edges are not mutual");

		assert( hfBase.face == hfNext.face && hfBase.face == hfNextNext.face && "Half edges are not in the same face" );
		assert( m_faceData[hfBase.face] ==i || m_faceData[hfBase.face] == hfBase.next || m_faceData[hfBase.face] == hfNext.next && "Half edges are not in the same face" );
	
		bool isVertPoint2HE(false);
		vvert_iterator itTemp;
		itTemp.m_cur = itTemp.m_end = &m_heData[hfBase.twin];
		do 
		{
			if (itTemp.m_cur->twin == m_vertData[hfBase.vert])
				isVertPoint2HE = true;

			assert (m_heData[itTemp.m_cur->twin].vert == hfBase.vert && "1-ring half edges don't share same vertex");
		} while (advance_iterator(itTemp));

		assert (isVertPoint2HE && "Center vertex points to other bunch of half edges");	
	}

	return true;
}


/*========================================
			Deformation
*========================================*/

void EditMesh::get_regionOfInf_vert(vector<size_t>& _vertVecInfluence)
{
	_vertVecInfluence.clear();
	
	for (auto iter(selected_reg_influence_verts.begin()); iter != selected_reg_influence_verts.end(); ++iter)
		_vertVecInfluence.push_back(*iter);
}
void EditMesh::clear_regionOfInf_vert()
{
	selected_reg_influence_verts.clear();
}

void EditMesh::get_control_vert(size_t& _control_vert)
{
	_control_vert = control_vert;
}

void EditMesh::clear_control_vert()
{
	 control_vert = std::numeric_limits<size_t>::max();
}

void EditMesh::clear_displacementQuanta()
{
	DisplacementQuanta[0] = DisplacementQuanta[1] = DisplacementQuanta[2] = std::numeric_limits<size_t>::max();
}

/*========================================
			computing weights
*========================================*/

double EditMesh::get_weight( const vvert_iterator& it ) const {
	const half_edge *itCur = it.m_cur;
	const half_edge *itTwin = &m_heData[ itCur->twin ];

	double result = 0;

	Eigen::Vector3d a = this->get_vertex( itTwin->vert );
	Eigen::Vector3d b = this->get_vertex( itCur->vert );
	
	Eigen::Vector3d e0 = (b - a);
	double eLen = e0.norm();

	e0 /= eLen;

	assert( ( itCur->face != HOLE_INDEX || itTwin->face != HOLE_INDEX ) && "Invalid mesh: edge with no face on either side" );

	if( itCur->face != HOLE_INDEX ){
		Eigen::Vector3d c = this->get_vertex( this->prev( *itCur ).vert );
		Eigen::Vector3d e1 = (c - a).normalized();

		result += std::tan( 0.5 * std::acos( e0.dot(e1) ) );
	}

	if( itTwin->face != HOLE_INDEX ){
		Eigen::Vector3d c = this->get_vertex( this->prev( *itTwin ).vert );
		Eigen::Vector3d e1 = (c - a).normalized();

		result += std::tan( 0.5 * std::acos( e0.dot(e1) ) );
	}
	
	return fabs(result / eLen);
}

double EditMesh::get_weight_sum(std::size_t vertex)
{
	double sum(0);

	// Boundary handling
	if (this->isBoundaryVert(vertex))
	{
		//std::cerr << "error because of boundary point" << std::endl;
		return std::numeric_limits<double>::max();
	}

	// initial the iterator visit the one ring of "vertex"
	vvert_iterator it;
	if (!this->init_iterator(it, vertex))
		std::cerr << "error " << std::endl;

	// sum the 1-ring angles	
	do {
		sum += this->get_weight(it);
		//sum += 1;
	} while(advance_iterator(it));

	return sum;	
}

double EditMesh::compWeights(std::size_t vFrom, std::size_t vTo ){
	//from j to i
	half_edge* he = this->find_edge(vFrom, vTo );
	Eigen::Vector3d left_v  = this->get_vertex(m_heData[m_heData[he->next].next].vert);
	Eigen::Vector3d right_v = this->get_vertex(m_heData[m_heData[m_heData[he->twin].next].next].vert);
	
	//compute alpha
	Eigen::Vector3d e0 = this->get_vertex(vFrom) - left_v;
	Eigen::Vector3d e1 = this->get_vertex(vTo) - left_v;
	// edge lengths of the two triangles, used to determine the 2d representation
	double    e0_len = e0.norm();
	double    e1_len = e1.norm();
	// We use the dot product to get cos 
	double  alpha = static_cast<double>(std::acos(( e0.dot(e1) ) / (e0_len*e1_len)));

	//compute betha
	Eigen::Vector3d e2 = this->get_vertex(vFrom) - right_v;
	Eigen::Vector3d e3 = this->get_vertex(vTo) - right_v;
	// edge lengths of the two triangles, used to determine the 2d representation
	double    e2_len = e2.norm();
	double    e3_len = e3.norm();
	// We use the dot product to get cos 
	double  betha = static_cast<double>(std::acos(( e2.dot(e3) ) / (e2_len*e3_len)));
	
	return 0.5*((std::cos(alpha)/std::sin(alpha))+(std::cos(betha)/std::sin(betha)));
	vvert_iterator it;	
	
}
/*========================================
			global and local solver
*========================================*/
Eigen::Vector3d EditMesh::get_old_delta(std::size_t vertex)
{
	double sumW(0);
	Eigen::Vector3d wieghted_vj_old(0, 0, 0);

	// initial the iterator visit the one ring of "vertex"
	vvert_iterator it;
	if (!this->init_iterator(it, vertex))
		std::cerr << "error" << std::endl;

	// sum the 1-ring angles	
	do {
		double weight(this->get_weight(it));
		//weight = 1;

		wieghted_vj_old += weight * m_vertices_old[it.m_cur->vert];

		sumW += weight;

	} while(advance_iterator(it));

	return (m_vertices_old[vertex]  - wieghted_vj_old / sumW);	
}

Eigen::Vector3d EditMesh::get_rotated_delta(std::size_t vertex)
{
	
	Eigen::Vector3d centerF(m_vertices_old[vertex]);
	Eigen::Vector3d centerT(m_vertices[vertex]);
	Eigen::Matrix3d crossCovariance = Eigen::Matrix3d::Zero();
	
	// compute centers
	size_t num(1);
	vvert_iterator it;
	if (!this->init_iterator(it, vertex))
		std::cerr << "error " << std::endl;	
	do {
		centerF += m_vertices_old[it.m_cur->vert];
		centerT   += m_vertices[it.m_cur->vert];

		++num;
	} while(advance_iterator(it));

	assert (num != 0);
	centerF /= num;
	centerT   /= num;

	// build cross variance matrix 3 * 3
	it.m_cur = it.m_end;
	Eigen::Vector3d tempFrom = m_vertices_old[vertex];
	Eigen::Vector3d tempTo = m_vertices[vertex];

	assert(tempFrom.rows() == 3 && tempFrom.cols() == 1 && tempTo.rows() == 3 && tempTo.cols() == 1 && 
		centerF.rows() == 3 && centerF.cols() == 1 && centerT.rows() == 3 && centerT.cols() == 1);
	crossCovariance += tempFrom * tempTo.transpose() - centerF * centerT.transpose();
	do {
		tempFrom = m_vertices_old[it.m_cur->vert];
		tempTo = m_vertices[it.m_cur->vert];

		assert(tempFrom.rows() == 3 && tempFrom.cols() == 1 && tempTo.rows() == 3 && tempTo.cols() == 1 && 
				centerF.rows() == 3 && centerF.cols() == 1 && centerT.rows() == 3 && centerT.cols() == 1);
		crossCovariance += tempFrom * tempTo.transpose() - centerF * centerT.transpose();
	} while(advance_iterator(it));

	crossCovariance /= num;

	// build LS matrix 4 * 4
	Eigen::Matrix3d tempA = crossCovariance - crossCovariance.transpose();
	Eigen::Vector3d delta(tempA(1, 2), tempA(2, 0), tempA(0, 1));	
	assert(delta.rows() == 3);
	Eigen::Matrix4d Arot(Eigen::Matrix4d::Zero());

	Arot.row(0) << crossCovariance.trace(), delta.transpose();
	Arot.col(0).tail(3) << delta; 
	Arot.block(1,1,3,3) << crossCovariance + crossCovariance.transpose() - crossCovariance.trace() * Eigen::Matrix3d::Identity();

	// get the smallest eigen vector
	Eigen::EigenSolver<Eigen::MatrixXd> es(Arot);
	double max(std::numeric_limits<double>::min());
	size_t position(0);
	for (size_t i(0); i < 4; ++i)
	{
		double lambda = es.eigenvalues()[i].real();
		if (max < lambda)
		{
			max = lambda;
			position = i;
		}
	}
	Eigen::Vector4d rot;
	rot[0] = es.eigenvectors().col(position)[0].real();
	rot[1] = es.eigenvectors().col(position)[1].real();
	rot[2] = es.eigenvectors().col(position)[2].real();
	rot[3] = es.eigenvectors().col(position)[3].real();
	rot.normalize();

	Eigen::Matrix3d rotMatrix;
	rotMatrix << rot[0] * rot[0] + rot[1] * rot[1] - rot[2] * rot[2] - rot[3] * rot[3], 
		2 * (rot[1] * rot[2] - rot[0] * rot[3]), 
		2 * (rot[1] * rot[3] + rot[0] * rot[2]),
		2 * (rot[1] * rot[2] + rot[0] * rot[3]),
		rot[0] * rot[0] + rot[2] * rot[2] - rot[1] * rot[1] - rot[3] * rot[3],
		2 * (rot[2] * rot[3] - rot[0] * rot[1]),
		2 * (rot[1] * rot[3] - rot[0] * rot[2]),
		2 * (rot[2] * rot[3] + rot[0] * rot[1]),
		rot[0] * rot[0] + rot[3] * rot[3] - rot[1] * rot[1] - rot[2] * rot[2];

	return rotMatrix * this->get_old_delta(vertex);
}

void EditMesh::popBmat(std::vector<size_t> ROI_verts, std::vector<bool> isBoundaryRegion,std::vector<size_t> ITL){
	// populate b = delta -> R[vi - (sum(w*vj_old))/sum(w)]
	for (size_t i(0), iEnd(ROI_verts.size()); i < iEnd; ++i)
	{
		if (isBoundaryRegion[i])
			continue;
	
		// compute residual 
		Eigen::Vector3d delta;
		if (m_vertices_old.size() != m_vertices.size())
		{
			m_vertices_old = m_vertices;
			//computing rotation including rotation
			delta = get_old_delta(ROI_verts[i]);
		}
		//the initial guess
		else if (IncludeRot)
			
			delta = get_rotated_delta(ROI_verts[i]);
		else
			delta = get_old_delta(ROI_verts[i]);
		
		size_t ind(ITL[i]);
		bX[ind] = delta[0];
		bY[ind] = delta[1];
		bZ[ind] = delta[2];
	}
}
std::vector<Eigen::Triplet<double>> EditMesh::popAmat(std::vector<size_t> ROI_verts, std::vector<bool> isBoundaryRegion,
											   std::vector<size_t> ITL, std::vector<size_t> ITR){
	// populate A

	std::vector<Eigen::Triplet<double>> A_elem;
	for (size_t i(0), iEnd(ROI_verts.size()); i < iEnd; ++i)
	{
		if (isBoundaryRegion[i])
			continue;
	
		// assign the diagonal entry
		assert (ITL[i] != HOLE_INDEX);
		size_t rowIndex(ITL[i]);
		A_elem.push_back(Eigen::Triplet<double>(rowIndex, rowIndex, 1) );
	
		// compute  delta and build the matrix A
		vvert_iterator it;	
		if (!this->init_iterator(it, ROI_verts[i]))
		{
			std::cerr << "error because of float point" << std::endl;
			break ;
		}
	

		double sum(this->get_weight_sum(ROI_verts[i]));
		double _sum(0);		// for debug
		do {
				double weight(this->get_weight(it));
				_sum += weight;


				size_t indexInRegionVec(ITR[it.m_cur->vert]);			
				assert(indexInRegionVec != HOLE_INDEX && indexInRegionVec < ROI_verts.size());
					
	
				if (it.m_cur->vert == control_vert)
				{
					Eigen::Vector3d temp(weight / sum * (m_vertices[control_vert] + DisplacementQuanta));

					bX[rowIndex] += temp[0];
					bY[rowIndex] += temp[1];
					bZ[rowIndex] += temp[2];
				}
				else if (isBoundaryRegion[indexInRegionVec])
				{
					Eigen::Vector3d temp(weight / sum * (m_vertices[it.m_cur->vert]));
	
					bX[rowIndex] += temp[0];
					bY[rowIndex] += temp[1];
					bZ[rowIndex] += temp[2];
				}
				else
				{
					size_t columnIndex(ITL[ITR[it.m_cur->vert]]);
					assert (columnIndex != rowIndex && columnIndex != HOLE_INDEX);

					A_elem.push_back(Eigen::Triplet<double>(rowIndex, columnIndex, - weight / sum));
				}
				
		} while(advance_iterator(it));
	}
	return A_elem;
}
void EditMesh::move_control_vert()
{
	//doing the necessary chhecks
	if (DisplacementQuanta[0] == std::numeric_limits<size_t>::max() &&
		 DisplacementQuanta[1] == std::numeric_limits<size_t>::max() &&
		  DisplacementQuanta[2] == std::numeric_limits<size_t>::max())
	{
		return;
	}

	if (DisplacementQuanta.norm() < 0.01)
		return;

	//ensure displacements are defined before doing anything
	assert( DisplacementQuanta[0] != std::numeric_limits<size_t>::max() && "error in computing displacement!" );
	assert( DisplacementQuanta[1] != std::numeric_limits<size_t>::max() && "error in computing displacement!" );
	assert( DisplacementQuanta[2] != std::numeric_limits<size_t>::max() && "error in computing displacement!" );


	//find the inner and the boundry of the region of influence
	// set the flag if a vert is inside the influence region
	std::vector<size_t> influence_region_verts;
	//contains vertices in the region of influence
	std::vector<bool> isInnerRegion(m_vertData.size(), false);	
	//containts vertices on the boundry of the region of influence 
	// plus the control vert
	std::vector<bool> isBoundaryRegion(selected_reg_influence_verts.size(), false);

	std::vector<size_t> indexToRegion(m_vertData.size(), HOLE_INDEX);
	//storing interior verts of the region of influence 
	std::vector<size_t> indexToLaplace(selected_reg_influence_verts.size(), HOLE_INDEX);

	for (auto iter(selected_reg_influence_verts.begin()); iter != selected_reg_influence_verts.end(); ++iter)
		influence_region_verts.push_back(*iter);

	// determine the verts that belong to the region of influence
	for (size_t i(0), iEnd(influence_region_verts.size()); i < iEnd; ++i)
	{
		isInnerRegion[influence_region_verts[i]] = true;
		indexToRegion[influence_region_verts[i]] = i;
	}

	// find the verts that are the boundary of the influence region
	for (size_t i(0), iEnd(influence_region_verts.size()); i < iEnd; ++i)
	{
		// treat the control vertex as a boudary vert of the region of influence
		if (influence_region_verts[i] == control_vert)
		{
			isBoundaryRegion[i] = true;
			continue;
		}
		
		// the "actual" boundary of the influence region
		vvert_iterator it;
		if (!this->init_iterator(it, influence_region_verts[i]))
		{
			std::cerr << "error" << std::endl;
			return;
		}

		do {
			if (!isInnerRegion[it.m_cur->vert])
			{
				isBoundaryRegion[i] = true;
				break;
			}
		} while(advance_iterator(it));
	}


	////////////////////////////////////////////////////////////////////////////
	//Main algorithm

	// set the index of interior verts for building Laplacian
	size_t index(0);
	for (size_t i(0), iEnd(influence_region_verts.size()); i < iEnd; ++i)
	{
		if (!isBoundaryRegion[i])
		{
			indexToLaplace[i] = index++;
		}
	}
	size_t InnerRegion_verts_n(index);

	//if only control vert is selcted
	if (InnerRegion_verts_n == 0)	// in this case we only move the control vertex as above
	{
		m_vertices[control_vert] += DisplacementQuanta;
		this->edit_count_updater();
		return;
	}

	//// slove the laplacian equation
	bX =Eigen::VectorXd::Zero(InnerRegion_verts_n);
	bY =Eigen::VectorXd::Zero(InnerRegion_verts_n);
	bZ =Eigen::VectorXd::Zero(InnerRegion_verts_n);
	Eigen::VectorXd Px(Eigen::VectorXd::Zero(InnerRegion_verts_n));
	Eigen::VectorXd Py(Eigen::VectorXd::Zero(InnerRegion_verts_n));
	Eigen::VectorXd Pz(Eigen::VectorXd::Zero(InnerRegion_verts_n));
	Eigen::SparseMatrix<double> A(InnerRegion_verts_n, InnerRegion_verts_n);

	// populate B
	popBmat(influence_region_verts,  isBoundaryRegion,indexToLaplace);

	// populate A
	std::vector<Eigen::Triplet<double>> A_elem = 
		popAmat(influence_region_verts, isBoundaryRegion, indexToLaplace, indexToRegion);

	A.setFromTriplets(A_elem.begin(), A_elem.end());
	
	// solve A^T * A * x = A^T * b using LLT
	bX = A.transpose() * bX;
	bY = A.transpose() * bY;
	bZ = A.transpose() * bZ;
	A = A.transpose() * A;
	Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(A);
	Px = solver.solve(bX);
	Py = solver.solve(bY);
	Pz = solver.solve(bZ);
	

	////m_vertices_old[control_vert] = m_vertices[control_vert];
	//////------------- The main algorithm 
	////m_vertices_old = m_vertices;
	////set_userConstraints(_selected_reg_influence_verts);
	////cout<<"done"<<endl;
	//////compAmat
	////Eigen::SparseMatrix<double> A = compAmat();
	//////comp  B_x, B_y, B_z
	////compBmat();

	////// solve using your prefered solver, we recommend LDLT or LLT
	//////new vert positions 
	//////A = A.transpose() * A;
	////Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(A.transpose() * A);
	////Eigen::VectorXd p_x,p_y,p_z;
	////p_x = solver.solve(A.transpose()*B_x);
	////p_y = solver.solve(A.transpose()*B_y);
	////p_z = solver.solve(A.transpose()*B_z);

	//////define m_vertices_new : the location of new vertices
	////for (size_t v=0;v<get_vert_number();v++){
	////	m_vertices_init.emplace_back( p_x[v], p_y[v], p_z[v] );
	////	//cout<<" vert : "<<v<<endl;
	////	//cout<<m_vertices[v].x()<<" "<<m_vertices_init[v].x()<<endl;
	////	//cout<<m_vertices[v].y()<<" "<<m_vertices_init[v].y()<<endl;
	////	//cout<<m_vertices[v].z()<<" "<<m_vertices_init[v].z()<<endl;
	////}
	//////do 3 iterations
	////for(size_t it=0;it<1;it++){
	////	cout<<"itteration "<<it<<endl;
	////	oneStep_iteration();
	////	m_vertices_init = m_vertices_new;
	////	m_vertices_new.clear();
	////	//cout<<"m_vertices_init "<<m_vertices_init.size()<<endl;
	////	//solve for rotation
	////
	////	/*cout<<"coordinate of original vs new vert positions are: "<<endl;
	////	for (size_t v=0;v<get_vert_number();v++){
	////		cout<<" vert : "<<v<<endl;
	////		cout<<m_vertices[v].x()<<" "<<m_vertices_init[v].x()<<endl;
	////		cout<<m_vertices[v].y()<<" "<<m_vertices_init[v].y()<<endl;
	////		cout<<m_vertices[v].z()<<" "<<m_vertices_init[v].z()<<endl;
	////	}
	////}

	////-------------------------------------------------------------------------


	//// update positions of vertices
	////for (size_t i(0), iEnd(_selected_reg_influence_verts.size()); i < iEnd; ++i)
	////{
	////	//m_vertices = m_vertices_init;
	////	if (!isBoundaryRegion[i])
	////		m_vertices[_selected_reg_influence_verts[i]] = Eigen::Vector3d(m_vertices_init[i].x(), m_vertices_init[i].y(), m_vertices_init[i].z());
	////}
	////


	// update positions of vertices
	for (size_t i(0), iEnd(influence_region_verts.size()); i < iEnd; ++i)
	{
		if (!isBoundaryRegion[i])
			m_vertices[influence_region_verts[i]] = Eigen::Vector3d(Px[indexToLaplace[i]], Py[indexToLaplace[i]], Pz[indexToLaplace[i]]);
	}
	
	// move the control vertex
	m_vertices[control_vert] += DisplacementQuanta;

	this->edit_count_updater();
}

Eigen::SparseMatrix<double> EditMesh::compAmat()
{
	

	//get number of rows for A,B, and X matrix
	cout<<" number of rows for A matrix: "<<get_rows()<<endl;
	//columns of A are euql to the number of verts
	Eigen::SparseMatrix<double> A (get_rows(), get_vert_number());
	std::vector<Eigen::Triplet<double>> A_elem;
	size_t i =0, iOld=0;

	for (size_t v=0;v<get_vert_number();v++)
	{
		std::vector<size_t> vNeighbor = find_borderVerts(v);
		
		for (size_t Vneigbor_number = 0;i<iOld+vNeighbor.size();i++,Vneigbor_number++)
		{
			double weigt = compWeights(vNeighbor[Vneigbor_number], v );
			double sum   = get_weight_sum (v);
			double w	 = weigt/sum;
			A_elem.push_back(Eigen::Triplet<double>(i,v, w) );
			A_elem.push_back(Eigen::Triplet<double>(i,vNeighbor[Vneigbor_number], -w) );
		}
		iOld += vNeighbor.size();
	}
	//set the soft constraints
	for(size_t k=0;i<A.rows();i++,k++)
		A_elem.push_back(Eigen::Triplet<double>(i,userConstraints[k]._vertI, softConst) );
	A.setFromTriplets(A_elem.begin(), A_elem.end());
	/*cout<<"A matrix: "<<endl;
	cout<<A<<endl;*/
	return A;

}

void EditMesh::compBmat(){
	
	//get number of rows for A,B, and X matrix
	//columns of B are euql to the number of verts
	 B_x = B_y = B_z = Eigen::VectorXd::Zero(get_rows());
	size_t i =0, iOld=0;

	for (size_t v=0;v<get_vert_number();v++)
	{
		std::vector<size_t> vNeighbor = find_borderVerts(v);
		
		for (size_t Vneigbor_number = 0;i<iOld+vNeighbor.size();i++,Vneigbor_number++)
		{
			double weigt = compWeights(vNeighbor[Vneigbor_number], v );
			double sum   = get_weight_sum (v);
			double w	 = weigt/sum;

			B_x[i] = w* (m_vertices[v].x()-m_vertices[vNeighbor[Vneigbor_number]].x());
			B_y[i] = w* (m_vertices[v].y()-m_vertices[vNeighbor[Vneigbor_number]].y());
			B_z[i] = w* (m_vertices[v].z()-m_vertices[vNeighbor[Vneigbor_number]].z());
		}
		iOld += vNeighbor.size();
	}
	//set the soft constraints
	for(size_t k=0;i<B_x.size();i++,k++){
		B_x[i] = softConst*userConstraints[k]._vertPos.x() ;
		B_y[i] = softConst*userConstraints[k]._vertPos.y() ;
		B_z[i] = softConst*userConstraints[k]._vertPos.z() ;
	}
	/*cout<<"B matrix: "<<endl;
	cout<<B_x<<endl;*/

}
size_t EditMesh::get_rows()
{
	//get number of rows for A,B, and X matrix
	//rows are euql to the number of neighbours of verts + user defined constraints
	size_t num_row =0;
	for (size_t v=0;v<get_vert_number();v++){
		std::vector<size_t> vNeighbor = find_borderVerts(v);
		num_row += vNeighbor.size();
	}
	num_row += userConstraints.size();
	return num_row;
}
//void EditMesh::set_userConstraints(std::vector<size_t> _selected_reg_influence_verts)
//{
//	userConstraints.clear();
//
//	for (size_t i(0), iEnd(_selected_reg_influence_verts.size()); i < iEnd; ++i)
//	{
//		if (_selected_reg_influence_verts[i] == control_vert)
//			continue;
//		
//		if (isBoundaryRegion[i])
//			userConstraints.push_back(userConstraint(_selected_reg_influence_verts[i], m_vertices[_selected_reg_influence_verts[i]]));
//	}
//	//for (size_t i(0), iEnd(_selected_reg_influence_verts.size()); i < iEnd; ++i)
//	//	{
//	//		if (isBoundaryRegion[i])
//	//				
//	userConstraints.push_back(userConstraint(control_vert,m_vertices[control_vert]));
//	//set the soft constraint
//	softConst = 1;
//}


//double EditMesh::get_weight_sum(std::size_t vertex)
//{
//	double sum(0);
//
//	// Boundary handling
//	if (this->isBoundaryVert(vertex))
//	{
//		//std::cerr << "error because of boundary point" << std::endl;
//		return std::numeric_limits<double>::max();
//	}
//
//	// initial the iterator visit the one ring of "vertex"
//	vvert_iterator it;
//	if (!this->init_iterator(it, vertex))
//		std::cerr << "error because of float point" << std::endl;
//
//	// sum the 1-ring angles	
//	do {
//		sum += this->compWeights(it.m_cur->vert, vertex );
//		//sum += 1;
//	} while(advance_iterator(it));
//
//	return sum;	
//}
Eigen::MatrixXd EditMesh::get_Arot(std::size_t v,std::size_t vj){
	//old vert position 3x9 matrix
	Eigen::MatrixXd Arot = Eigen::MatrixXd::Zero(3,9);
	size_t j=0, j_old =0;
	
	for(size_t i=0;i<Arot.rows();i++,j++){
		/*Arot(i,j) = (this->get_vertex(v).x()-this->get_vertex(vj).x());
		j++;
		Arot(i,j) = (this->get_vertex(v).y()-this->get_vertex(vj).y());
		j++;
		Arot(i,j) = (this->get_vertex(v).z()-this->get_vertex(vj).z());*/

		Arot(i,j) = (m_vertices_old[v].x()-m_vertices_old[vj].x());
		j++;
		Arot(i,j) = (m_vertices_old[v].y()-m_vertices_old[vj].y());
		j++;
		Arot(i,j) = (m_vertices_old[v].z()-m_vertices_old[vj].z());
	}
	/*cout<<"A matrix for computing rotation: "<<endl;
	cout<<Arot<<endl;*/
	return Arot;

}
Eigen::MatrixXd EditMesh::get_R(std::size_t v){

	//get neighboring verts (CW)
	vector<size_t> border_verts = find_borderVerts(v);
	

	//compute sum(LSrot_vectArot_i'*Arot_i) for the vert and its neigbors 
	vector<Eigen::MatrixXd>LSrot_vect;
	vector<Eigen::MatrixXd>RSrot_vect, RSrot2_vect;
	Eigen::VectorXd vert_new = Eigen::VectorXd::Zero(3,1);
	Eigen::MatrixXd sum_LSrot, sum_RSrot;
	sum_LSrot=  Eigen::MatrixXd::Zero(9,9), sum_RSrot=  Eigen::MatrixXd::Zero(9,1);
	for(size_t i=0;i<border_verts.size();i++){

		double weigt = compWeights(border_verts[i], v );
		double sum   = get_weight_sum (v);
		double w	 = weigt/sum;
		
		//compute the lhs
		LSrot_vect.push_back(	w*(this->get_Arot(v,border_verts[i]).transpose())	*
			(this->get_Arot(v,border_verts[i]))	);
		//compute the rhs
		RSrot_vect.push_back(	w*(this->get_Arot(v,border_verts[i]).transpose())	*
			(m_vertices_init[v] - m_vertices_init[border_verts[i]])	);
		
		sum_LSrot += LSrot_vect[i];
		sum_RSrot += RSrot_vect[i];
	}

	Eigen::VectorXd R_approx = (sum_LSrot).ldlt().solve(sum_RSrot);
	Eigen::MatrixXd Rmat_approx(3,3) ;
	Rmat_approx <<R_approx(0), R_approx(1), R_approx(2),
				  R_approx(3), R_approx(4), R_approx(5),
				  R_approx(6), R_approx(7), R_approx(8);
	//SVD portion
	//Rmat_approx = USV*
	//The solution for R is the product UV∗
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Rmat_approx, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXd R = svd.matrixU()*svd.matrixV();
	//cout<<"Rotation matrix R is: "<<endl;
	//cout<<R<<endl;
	return R;
}
void EditMesh::oneStep_iteration(){
	//----------compute the lhs-------------------//
	Eigen::SparseMatrix<double> A = compAmat();
	//----------compute the rhs-------------------//
	//get number of rows for A,B, and X matrix
	//columns of B are euql to the number of verts
	size_t i =0, iOld=0;
	Eigen::VectorXd Rs_x, Rs_y, Rs_z;
	Rs_x = Rs_y = Rs_z= Eigen::VectorXd::Zero(get_rows());
;
	for (size_t v=0;v<get_vert_number();v++)
	{
		std::vector<size_t> vNeighbor = find_borderVerts(v);
		

		for (size_t Vneigbor_number = 0;i<iOld+vNeighbor.size();i++,Vneigbor_number++)
		{
			this->get_R(v);
			double weigt = compWeights(vNeighbor[Vneigbor_number], v );
			double sum   = get_weight_sum (v);
			double w	 = weigt/sum;

			Rs_x[i] = w/2 * (this->get_R(v).row(0) + this->get_R(vNeighbor[Vneigbor_number]).row(0)) *
				(m_vertices[v] - m_vertices[vNeighbor[Vneigbor_number]]);
			Rs_y[i] = w/2 * (this->get_R(v).row(1) + this->get_R(vNeighbor[Vneigbor_number]).row(1)) *
				(m_vertices[v]-m_vertices[vNeighbor[Vneigbor_number]]);
			Rs_z[i] = w/2 * (this->get_R(v).row(2) + this->get_R(vNeighbor[Vneigbor_number]).row(2)) *
				(m_vertices[v]-m_vertices[vNeighbor[Vneigbor_number]]);
			
		}
		iOld += vNeighbor.size();
	}
	//set the soft constraints //this portion doesn't change -> there is no rotation
	for(size_t k=0;i<Rs_x.size();i++,k++){
		Rs_x[i] = softConst*userConstraints[k]._vertPos.x() ;
		Rs_y[i] = softConst*userConstraints[k]._vertPos.y() ;
		Rs_z[i] = softConst*userConstraints[k]._vertPos.z() ;
	}
	
	// solve using your prefered solver, we recommend LDLT or LLT
	//to solvde for new vert positions 
	Eigen::VectorXd v_x, v_y, v_z;
	Rs_x = A.transpose() * Rs_x;
	Rs_y = A.transpose() * Rs_y;
	Rs_z = A.transpose() * Rs_z;
	A = A.transpose() * A;
	Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > solver(A);
	v_x = solver.solve(Rs_x);
	v_y = solver.solve(Rs_y);
	v_z = solver.solve(Rs_z);

	//define m_vertices_new : the location of new vertices
	for (size_t v=0;v<get_vert_number();v++){
		m_vertices_new.emplace_back( v_x[v], v_y[v], v_z[v] );
	}
}


/*========================================
			Remeshing
*========================================*/
void EditMesh::writeTrianglesAndVertsToFile()
{
	std::ofstream myfile;
	myfile.open("generated_data_30/xyzPositions.txt");
	for (int i = 0; i < vertsInput.size(); i++)
		myfile << vertsInput[i] << "\n";
	myfile.close();


	myfile.open("generated_data_30/triangleVerts.txt");
	for (int i = 0; i < triangles.size(); i++)
		myfile << triangles[i] << "\n";
	myfile.close();

}

void EditMesh::readTrianglesAndVertsFromFileAndBuildMesh56()
{
	std::string line;
	
	std::ifstream myfile("generated_data_56/xyzPositions.txt");
	//std::ifstream myfile(xyzPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num = std::stod(line, &numS);
			//std::cout << num << '\n';
			vertsInput.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file xyzPositions";

	myfile = std::ifstream("generated_data_56/triangleVerts.txt");
	//myfile = std::ifstream(trianglesPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num =  std::stod(line, &numS);
			//std::cout << num << '\n';
			triangles.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file triangleVerts";

	this->init(vertsInput, triangles);
	edit_count++;
}

void EditMesh::readTrianglesAndVertsFromFileAndBuildMesh113()
{
	std::string line;
	
	std::ifstream myfile("generated_data_113/xyzPositions.txt");
	//std::ifstream myfile(xyzPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num = std::stod(line, &numS);
			//std::cout << num << '\n';
			vertsInput.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file xyzPositions";

	myfile = std::ifstream("generated_data_113/triangleVerts.txt");
	//myfile = std::ifstream(trianglesPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num =  std::stod(line, &numS);
			//std::cout << num << '\n';
			triangles.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file triangleVerts";

	this->init(vertsInput, triangles);
	edit_count++;
}

void EditMesh::readTrianglesAndVertsFromFileAndBuildMesh30()
{
	std::string line;
	
	std::ifstream myfile("generated_data_30/xyzPositions.txt");
	//std::ifstream myfile(xyzPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num = std::stod(line, &numS);
			//std::cout << num << '\n';
			vertsInput.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file xyzPositions";

	myfile = std::ifstream("generated_data_30/triangleVerts.txt");
	//myfile = std::ifstream(trianglesPath);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			std::string::size_type numS;
			double num =  std::stod(line, &numS);
			//std::cout << num << '\n';
			triangles.push_back(num);
		}
		myfile.close();
	}
	else std::cout << "Unable to open file triangleVerts";

	this->init(vertsInput, triangles);
	edit_count++;
}

void EditMesh::rm_flip_edge(){
	size_t num_edge_flips =0;
	for (int i = 0; i<m_heData.size();i++)
	{
		//half-edge to flip
		 half_edge &he = m_heData[i];
		 half_edge &twin = m_heData[ he.twin ];

		 //compute the old diagonal before edge flip
		 Eigen::Vector3d vFrom_old = this->get_vertex( he.vert );
		 Eigen::Vector3d vTo_old   = this->get_vertex( twin.vert );
		 if( he.face == HOLE_INDEX || twin.face == HOLE_INDEX )
			 continue;
		 if( this->isBoundaryEdge(twin.twin))
			 continue;

		 double eO_len = (vTo_old -vFrom_old).norm();
		 
		 //compute the new diagonal after edge flip
		 size_t he_New   = m_heData[ he.next ].next;
		 size_t twin_New = m_heData[twin.next].next;
		 Eigen::Vector3d vFrom_new = this->get_vertex(m_heData[he_New].vert );
		 Eigen::Vector3d vTo_new   = this->get_vertex(m_heData[twin_New].vert );
		 double eNew_len = (vTo_new -vFrom_new).norm();

		 //flip if old diagonal longer than the new one
		 if (eO_len > eNew_len)
		 {
			 num_edge_flips++;
			 //test normal error//
			 //vert notmals from original mesh
			 
			 // prep: gather half edge indices in
			// the order they should be after flip
			std::size_t he_tri[3],twin_tri[3];
			he_tri[0]   = he.next;
			twin_tri[0] = twin.next;
			he_tri[1]   = twin.twin;
			twin_tri[1] = he.twin;
			he_tri[2]   = m_heData[ twin_tri[0] ].next;
			twin_tri[2] = m_heData[ he_tri[0] ].next;

			// step 2: set the he's vert to new originating vert
			size_t Vhe_tri1   = m_heData[ twin_tri[2] ].vert;
			size_t Vtwin_tri1 =  m_heData[ he_tri[2] ].vert;
			
			// step 3: compute the vet normals from the original mesh
			Eigen::Vector3d v_tri_norm[3],vTwin_tri_norm[3];
			Eigen::Vector3d v_tri[3],vTwin_tri[3];

			/*v_tri_norm[0]	  =  this-> get_vnormal(m_heData[he_tri[0]].vert);
			v_tri_norm[1]	  =  this-> get_vnormal(Vhe_tri1);
			v_tri_norm[2]	  =  this-> get_vnormal(m_heData[he_tri[2]].vert);
			vTwin_tri_norm[0] =  this-> get_vnormal(m_heData[twin_tri[0]].vert);
			vTwin_tri_norm[1] =  this-> get_vnormal(Vtwin_tri1);
			vTwin_tri_norm[2] =  this-> get_vnormal(m_heData[twin_tri[2]].vert);*/

			v_tri_norm[0]	  =  oldVertices_normals[m_heData[he_tri[0]].vert];
			v_tri_norm[1]	  =  oldVertices_normals[m_heData[he_tri[1]].vert];
			v_tri_norm[2]	  =  oldVertices_normals[m_heData[he_tri[2]].vert];
			vTwin_tri_norm[0] =  oldVertices_normals[m_heData[twin_tri[0]].vert];
			vTwin_tri_norm[1] =  oldVertices_normals[Vtwin_tri1];
			vTwin_tri_norm[2] =  oldVertices_normals[m_heData[twin_tri[2]].vert];

			// step 4: compute the current face normals 
			v_tri[0]	 = this->get_vertex( m_heData[he_tri[0]].vert );
			v_tri[1]	 = this->get_vertex( Vhe_tri1 );
			v_tri[2]	 = this->get_vertex( m_heData[he_tri[2]].vert );
			vTwin_tri[0] = this->get_vertex( m_heData[twin_tri[0]].vert );
			vTwin_tri[1] = this->get_vertex( Vtwin_tri1 );
			vTwin_tri[2] = this->get_vertex( m_heData[twin_tri[2]].vert );
			Eigen::Vector3d fnormal_tri	   = ((v_tri[1] - v_tri[0]).cross(v_tri[2] - v_tri[0])).normalized();
			Eigen::Vector3d fnormalTwin_tri = ((vTwin_tri[1] - vTwin_tri[0]).cross(vTwin_tri[2] - vTwin_tri[0])).normalized();

			double errorLimit = 1.0 - std::cos(PI/9); //10 degrees
			bool errorFlag = true;
			// step 5: compute the error (dot product between old verts and the current faces
			for(int i=0;i<3;i++){
				double error_tri	 =  1.0 - v_tri_norm[i].dot(fnormal_tri);
				double errorTwin_tri =  1.0 - vTwin_tri_norm[i].dot(fnormalTwin_tri);
				if(error_tri >errorLimit || errorTwin_tri>errorLimit){
					errorFlag = false;
					break;
				}
			}
			if (errorFlag == false)
				continue;

			this->flip_edge(he);
			
		 }
	
	}
	cout<<num_edge_flips<<" out of "<<m_heData.size()<<"edges  are flipped"<<endl;
}

void EditMesh::refine_mesh(){
	size_t j =0;
	for (int i = 0; i<m_heData.size();i++)
	{
		//half-edge to refine
		 half_edge &he = m_heData[i];
		 half_edge &twin = m_heData[ he.twin ];
		 //compute the old diagonal before edge flip
		 Eigen::Vector3d vFrom_old = this->get_vertex( he.vert );
		 Eigen::Vector3d vTo_old   = this->get_vertex( twin.vert );

		 if( he.face == HOLE_INDEX || twin.face == HOLE_INDEX )
			 continue;

		 double e_len = (vTo_old -vFrom_old).norm();
		 //cout<<"e_len"<<e_len<<endl;
		 double elenlimit = 0.13; //this could change, just some arbitrary value 
		
		 if(e_len>elenlimit)
		 {
			 
			 //location of new vert
			 size_t newVertI = m_vertices.size();
			 Eigen::Vector3d newVertPosition(vFrom_old+((vTo_old -vFrom_old)/2.0));

			 ///cout<<"newVertPosition: "<<newVertPosition.x()<<" "<<newVertPosition.y()<<" "<<newVertPosition.z()<<" "<<endl;
			
			 //step 1: prep: gather half edge indices in
			// the order they should be after adding the vert
			std::size_t  he_tri_lu[3], he_tri_ru[3], he_tri_ld[3], he_tri_rd[3];
			half_edge heNew[6];
			//old hes
			he_tri_lu[0] = twin.twin;
			he_tri_lu[1] = m_heData.size();
			he_tri_lu[2] = m_heData[he.next].next;
			he_tri_ld[0] = he.twin;
			he_tri_ld[1] = twin.next;
			he_tri_ld[2] = m_heData.size()+1;
			he_tri_ru[0] = m_heData.size()+2;
			he_tri_ru[1] = he.next;
			he_tri_ru[2] = m_heData.size()+3;
			he_tri_rd[0] = m_heData.size()+4;
			he_tri_rd[1] = m_heData.size()+5;
			he_tri_rd[2] = m_heData[ twin.next ].next;

			// Update the face index of the right two triangles.
			size_t faceIndex_ru = m_faceData.size();
			m_faceData.push_back(he_tri_ru[1]);
			size_t faceIndex_rd = m_faceData.size();
			m_faceData.push_back(he_tri_rd[2]);
			//make sure the faces on the left side map to one of the hes on the left that will 
			//stll remain on the boundry of the face after the flip
			m_faceData[m_heData[he_tri_lu[0]].face] = he_tri_lu[0];
			m_faceData[m_heData[he_tri_ld[0]].face] = he_tri_ld[0];

			//add the new vert
			m_vertices.push_back(newVertPosition);
			// A mapping from vertex index to an arbitrary half-edge originating from this vertex.
			m_vertData.push_back(he_tri_lu[1]);
			

			//step 2: create 6 new hes:
			//(1)lu[1] (2)ld[2] (3)ru[0] (4)ru[2] (5)rd[0] (6)rd[1]

			//he_tri_lu[1]:need to set the twin
			detail::init( heNew[0], he_tri_lu[2],he_tri_ru[2], newVertI,  m_heData[he_tri_lu[0]].face );
			//he_tri_ld[2]:need to set the twin
			detail::init( heNew[1], he_tri_ld[0], he_tri_rd[1], m_heData[ he_tri_rd[2]].vert,  m_heData[he_tri_ld[0]].face );
			//he_tri_ru[0]:need to set the twin
			detail::init( heNew[2], he_tri_ru[1], he_tri_rd[0], newVertI,  faceIndex_ru);
			//he_tri_ru[2]:need to set the twin, next
			detail::init( heNew[3], he_tri_ru[0], he_tri_lu[1], m_heData[he_tri_lu[2]].vert,  faceIndex_ru );
			//he_tri_rd[0]:need to set the twin, next
			detail::init( heNew[4],he_tri_rd[1], he_tri_ru[0], m_heData[he_tri_ru[1]].vert,  faceIndex_rd );
			//he_tri_rd[1]:need to set the twin
			detail::init( heNew[5], he_tri_rd[2], he_tri_ld[2], newVertI,  faceIndex_rd );
			
			//step 3:update the necessary attributes of the old hes
			//(1) he_tri_lu[0]
			m_heData[he_tri_lu[0]].next = he_tri_lu[1];
			//(2) he_tri_lu[2] -> it's fine as it is
			//(3) he_tri_ld[0]
			m_heData[he_tri_ld[0]].vert = newVertI;
			//(4) he_tri_ld[1]
			m_heData[he_tri_ld[1]].next = he_tri_ld[2];
			//(5) he_tri_ru[1]
			m_heData[he_tri_ru[1]].next = he_tri_ru[2];
			m_heData[he_tri_ru[1]].face = faceIndex_ru;
			//(6) he_tri_ru[1]
			m_heData[he_tri_rd[2]].next = he_tri_rd[0];
			m_heData[he_tri_rd[2]].face = faceIndex_rd;

			//step 4: push back the new collected data 
			m_heData.push_back(heNew[0]);
			m_heData.push_back(heNew[1]);
			m_heData.push_back(heNew[2]);
			m_heData.push_back(heNew[3]);
			m_heData.push_back(heNew[4]);
			m_heData.push_back(heNew[5]);
			
			j++;
			
		 } 
	}
	cout<<j<<" vertices were added"<<endl;
	newVertsNum = j;
}

//Move vertices ON surface to improve
//sizing/quality
bool EditMesh::get_smoothnessError(size_t vertI){
	vvert_iterator it;
	double errorLimit = 1.0 - std::cos(PI/18); //10 degrees
	if (!this->init_iterator(it, vertI))
		std::cerr << "error because of float point" << std::endl;

	// 1-ring hes originating from this vert	
	do {
		 Eigen::Vector3d fnormal = get_fnormal(it.m_cur->face);

		 double error_tri	 =  1.0 - this-> get_vnormal(vertI).dot(fnormal);
		 if(error_tri >errorLimit)
			 return false;
	} while(advance_iterator(it));
	return true;
}

void EditMesh::smooth_mesh()
{
	size_t numVertsMoved =0;
	for(size_t i=0;i<m_vertData.size();i++)
	{
		//check after you simplify
		if(m_vertData[i] == HOLE_INDEX)
			continue;

		//check to ensure the vert we are moving is not a boundary vert
		if(this->isBoundaryVert(i))
			continue;

		std::vector<size_t> umbrellaVerts = find_borderVerts(i);
		//assuming equal weights on each of the border verts
		Eigen::Vector3d newVertPos, oldVertPos;
		Eigen::Vector3d sum;
		sum.setZero();
		for(int it=0;it<umbrellaVerts.size();it++)
			sum += this->get_vertex(umbrellaVerts[it]);
		
		newVertPos = sum/umbrellaVerts.size();
		oldVertPos = this->get_vertex(i);

		//check the error metrics
		//check the original vert normal with its umbrella face normals
		//if(!get_smoothnessError(i))
			//continue;
		
		//compute normal-based error
		//get the old and the new verts normals
		Eigen::Vector3d vOldNormal, vNewNormal;
		vOldNormal = oldVertices_normals[i];
		m_vertices[i] = newVertPos;
		vNewNormal = this-> get_vnormal(i);
		//compute the error
		double errorLimit = 1.0 - std::cos(PI/18); //10 degrees
		double errorVert = 1.0 - vOldNormal.dot(vNewNormal);

		//cout<<"error limit is: "<<errorLimit<<endl;
		//cout<<"errorVert is "<<errorVert<<endl;
		if(errorVert>errorLimit){
			m_vertices[i] = oldVertPos;
			continue;
		}
		if(!get_smoothnessError(i)){
			m_vertices[i] = oldVertPos;
			continue;
		}
		
		numVertsMoved++;
	}
	cout<<"in total "<<numVertsMoved<<" out of "<<m_vertices.size()<<" vetices were moved."<<endl;
}
void EditMesh::re_mesh(){

	/*//simplify 500 times
	cout<<m_vertData.size()<<endl;
	for (int j =0;j<100;j++){
		collapse_one_edge(*this);
		//this->oneStep_simplify();
	}
	
	int vertSize =0;
	for (int k=0;k<m_vertData.size();k++){
		if(m_vertData[k] == HOLE_INDEX)
			continue;
		vertSize++;
	}
	cout<<vertSize<<endl;
	*/
	


	for(int i=0;i<m_vertData.size();i++)
		oldVertices_normals.push_back(this->get_vnormal(i));

	/*computeInitialErrors_advancedMethod();
	sortVerticesByError();
	numOfRemovedVert = 0;
	simplify(1000);*/
	
	for(size_t it=0;it<100;it++){

		//collapse_one_edge(*this);
		//add vertices
		this->refine_mesh();

		//after every refinement you need to update the oldvertnormals
		for (size_t k=newVertsNum; k>0;k--){
			oldVertices_normals.push_back(this->get_vnormal(m_vertices.size()-k));
		}
		/*for (size_t k=newVertsNum; k>0;k--){
			collapse_one_edge(*this);
		}*/
		
		//flip an edge
		this->rm_flip_edge();
		
		//move a vertex
		this->smooth_mesh();
		
		//flip an edge
		this->rm_flip_edge();
	}
	
	this->edit_count_updater();	// for mesh view
}

void EditMesh::loadMesh(){

	ImageDataLoader test;
	test.numImages = 32;
	//test.numImages = 7;

	//std::vector<Eigen::MatrixXd> images = test.readAllImages("lattice_test_images_data", test.numImages);
	//std::vector<Eigen::MatrixXd> images = test.readAllImages("test_images_data", test.numImages);
	std::vector<Eigen::MatrixXd> images = test.readAllImages("images_data_denoised", test.numImages);

	//for (int i = 0; i < images.size(); i++)
	//{
	//	cout << test.imagePath[i] << "\n" << images[i] << "\n";
	//}
	
	createImagesGrid(images);
	createGridCells(images);

	this->writeTrianglesAndVertsToFile();
	//this->readTrianglesAndVertsFromFileAndBuildMesh();

}

void EditMesh::load30Slices(){
	readTrianglesAndVertsFromFileAndBuildMesh30();
}
void EditMesh::load56Slices(){
	readTrianglesAndVertsFromFileAndBuildMesh56();
}
void EditMesh::load113Slices(){
	readTrianglesAndVertsFromFileAndBuildMesh113();
}
//-------------------------------------- Simplification - Yasmin
#pragma region simplification
void EditMesh::simplify(size_t numOfVertsToRemove)
{
	while (numOfVertsToRemove != numOfRemovedVert && !isSimplest())
	{
		size_t vi;

		vi = getLeastErrorVertex();

		if (!isVertexRemovalValid(vi)){

			moveVertextToEnd(vi);

			//printVerticesWithErrors();

			vi = getLeastErrorVertex();

		}
		else
		{
			//printf("\nValid. Using vi = [%d]\n", vi);

			distributeErrorToNeighbors(vi);

			//printVerticesWithErrors();

			deleteFacesAndTriangulate(vi);

			sortVerticesByError();

			numOfRemovedVert++;

		}
	}

	numOfRemovedVert = 0;

	edit_count++;
}

bool EditMesh::isSimplest()
{
	if (vertices.size() == 4)
		return true;
	if (vertices.at(vertices.size() - 5).error == HOLE_INDEX)
		return true;
	else
		return false;
}

void EditMesh::moveVertextToEnd(size_t vi)
{
	size_t insertlocation = 0;
	for (int i = 0; i < vertices.size(); i++)
	if (vertices[i].error != HOLE_INDEX)
	{
		insertlocation = i;
		break;
	}
	for (int i = 0; i < vertices.size(); i++)
	if (vi == vertices[i].vertexIndex)
	{
		vertices_error v = vertices[i];
		vertices.pop_back();
		vertices.insert(vertices.begin() + insertlocation, v);
	}
}

bool EditMesh::isVertexRemovalValid(size_t vi)
{
	std::vector<size_t> bv = getBorderVertices(vi);

	size_t ov = bv[0];

	half_edge* he;

	size_t connV1i; //prev vertex
	size_t connV2i; //next vertex

	for (int i = 0; i < bv.size(); i++)	{
		if (i == 0)
		{
			connV1i = i + 1;
			connV2i = bv.size() - 1;
		}
		else if (i == bv.size() - 1)
		{
			connV1i = i - 1;
			connV2i = 0;
		}
		else
		{
			connV1i = i - 1;
			connV2i = i + 1;
		}

		//printf("i = %d, size = %d, connV1i = %d, connV2i = %d\n", i, bv.size(), connV1i, connV2i);

		for (int j = 0; j < bv.size(); j++)
		{
			if (bv[j] != bv[i] && bv[j] != bv[connV1i] && bv[j] != bv[connV2i])
			{
				he = find_edge(bv[i], bv[j]);
				if (he != NULL)
				{
					//printf("ERROR: i = %d, j = %d\n", i, j);
					//printf("\nNON-MANIFOLD [CASE 1]. Using vi = %d.\n", vi);
					return false;

				}
			}
		}

	}

	// check second case
	for (int i = 0; i < bv.size(); i++)
	{
		ov = bv[i];

		half_edge* heBase = find_edge(ov, vi);
		half_edge heTwin = m_heData[heBase->twin];

		// We are going to delete the faces on either side of the chosen edge, so we need to delete 3 half_edges and patch up the twin links on the 4
		// bordering edges.
		std::size_t heBorder[4];
		heBorder[0] = m_heData[heBase->next].twin;
		heBorder[1] = m_heData[m_heData[heBase->next].next].twin;
		heBorder[2] = m_heData[m_heData[heTwin.next].next].twin;
		heBorder[3] = m_heData[heTwin.next].twin;

		// Check if we can actually collapse. This checks for a degree 3 vertex at the vertices not on the edge we are collapsing.
		if (m_heData[m_heData[m_heData[heBorder[1]].next].twin].next == heBorder[0])
		{
			//printf("\nNON-MANIFOLD [CASE 2]. Using vi = %d.\n", vi);
			return false;
		}

		if (m_heData[m_heData[m_heData[heBorder[2]].next].twin].next == heBorder[3])
		{
			//printf("\nNON-MANIFOLD [CASE 2]. Using vi = %d.\n", vi);
			return false;
		}
	}

	return true;

}

size_t EditMesh::getNumRemainingVerts()
{
	size_t remainingVerts = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].error != HOLE_INDEX)
			remainingVerts++;
	}

	return remainingVerts;
}

size_t EditMesh::getLeastErrorVertex()
{
	//return vertices[vertices.size()-3].vertexIndex;
	return vertices.back().vertexIndex;
}

void EditMesh::computeInitialErrors_simpleMethod()
{
	int verticesNum = m_vertData.size();

	float error;
	for (int i = 0; i < verticesNum; i++)
	{
		error = computeErrorAtVertex(i);

		vertices_error newVert;
		newVert.error = error;
		newVert.vertexIndex = i;
		vertices.push_back(newVert);
	}
}

void::EditMesh::computeInitialErrors_advancedMethod()
{
	//distance to the approximating plane of the neighboring vertices 

	int verticesNum = m_vertData.size();

	float error;
	for (int i = 0; i < verticesNum; i++)
	{
		error = 0;
		vertices_error newVert;
		newVert.error = error;
		newVert.vertexIndex = i;
		vertices.push_back(newVert);
	}

	std::vector<size_t> bv;
	Eigen::Vector3d averageLocation;
	Eigen::Vector3d sumLocation;
	Eigen::Vector3d normal;
	float distanceFromUmbrellaMid;

	for (int i = 0; i < vertices.size(); i++)
	{
		bv = getBorderVertices(vertices[i].vertexIndex);

		for (int j = 0; j < bv.size(); j++)
			sumLocation = m_vertices[vertices[bv[j]].vertexIndex];

		averageLocation = sumLocation / bv.size();

		normal = m_vertices[vertices[i].vertexIndex] - averageLocation;

		distanceFromUmbrellaMid = sqrt(pow(normal.x(), 2) + pow(normal.y(), 2) + pow(normal.z(), 2));

		vertices[i].error = distanceFromUmbrellaMid;
	}

	printVerticesWithErrors();

}

void EditMesh::deleteFacesAndTriangulate(size_t vi)
{
	// get border vertices
	std::vector<size_t> borderVertices = getBorderVertices(vi);

	deleteUmbrellaFacesY(vi);

	triangulateHole(vi, borderVertices);

	//m_vertData[vi] = HOLE_INDEX;
}

std::vector<size_t> EditMesh::getBorderVertices(size_t vertIndex)
{
	std::vector<size_t> bVertices;
	std::vector<half_edge> connectedHe = getConnectedHalfEdges(vertIndex);

	for (int i = 0; i < connectedHe.size(); i++)
		bVertices.push_back(connectedHe[i].vert);

	return bVertices;
}

std::vector<half_edge> EditMesh::getBorderHalfEdges(std::vector<size_t> borderingVertices)
{
	std::vector<half_edge> b;

	size_t ov;
	std::vector<size_t> v(borderingVertices.size() - 1);

	ov = borderingVertices[0];
	for (int i = 0; i < borderingVertices.size() - 1; i++)
	{
		v[i] = borderingVertices[i + 1];
	}

	std::vector<half_edge> c = getConnectedHalfEdges(ov);

	for (int i = 0; i < c.size(); i++)
	if (c[i].face == HOLE_INDEX) b.push_back(c[i]);

	for (int i = 0; i < v.size(); i++)
	{
		c = getConnectedHalfEdges(v[i]);
		for (int j = 0; j < c.size(); j++)
		if (c[j].face == HOLE_INDEX) b.push_back(c[j]);
	}

	//for (int i = 0; i < b.size(); i++)
	//printf("half edge [%d]: f:%d, n:%d, t:%d, v:%d\n", i, b[i].face, b[i].next, b[i].twin, b[i].vert);

	return b;
}

std::vector<size_t> EditMesh::getBorderHalfEdgesIndices(std::vector<size_t> borderingVertices)
{
	std::vector<size_t> b;

	size_t ov;
	std::vector<size_t> v(borderingVertices.size() - 1);

	ov = borderingVertices[0];
	for (int i = 0; i < borderingVertices.size() - 1; i++)
	{
		v[i] = borderingVertices[i + 1];
	}

	std::vector<size_t> c = getConnectedHalfEdgesIndices(ov);

	for (int i = 0; i < c.size(); i++)
	if (m_heData[c[i]].face == HOLE_INDEX) b.push_back(c[i]);

	for (int i = 0; i < v.size(); i++)
	{
		c = getConnectedHalfEdgesIndices(v[i]);
		for (int j = 0; j < c.size(); j++)
		if (m_heData[c[j]].face == HOLE_INDEX) b.push_back(c[j]);
	}

	//for (int i = 0; i < b.size(); i++)
	//	printf("half edge [%d]: f:%d, n:%d, t:%d, v:%d\n", i, m_heData[b[i]].face, m_heData[b[i]].next,
	//	m_heData[b[i]].twin, m_heData[b[i]].vert);

	return b;
}

void EditMesh::deleteUmbrellaFacesY(size_t vertIndex)
{
	int he;
	while (m_vertData[vertIndex] != HOLE_INDEX)
	{
		he = m_vertData[vertIndex];
		delete_face(m_heData[he].face);
	}
}

void EditMesh::printVerticesWithErrors()
{
	for (int i = 0; i< vertices.size(); i++)
		printf("[%d: %d, %1.2f]\n ", i, vertices.at(i).vertexIndex, vertices.at(i).error);
	printf("==================\n");
}

void EditMesh::triangulateHole(size_t vertIndex, std::vector<size_t> borderVertices)
{

	// get the indices of the border half edges
	std::vector<size_t> bIndex = getBorderHalfEdgesIndices(borderVertices);


	// get number of border edges
	size_t n = bIndex.size();
	// get border vertices (v) and originating vertex (ov)
	std::vector <size_t> v;
	size_t ov;
	for (int i = 0; i < bIndex.size(); i++)
	{
		if (i == 0)
			ov = m_heData[bIndex[i]].vert;
		else
			v.push_back(m_heData[bIndex[bIndex.size() - i]].vert);
	}

	//printf("\nov = %d", ov);
	//for (int i = 0; i < v.size(); i++)printf(", v = %d", v[i]); printf("\n");

	// get base indices for the new faces and half edges to be created
	size_t baseFaceIndex = m_faceData.size();
	size_t heBaseIndex = m_heData.size();

	//printf("baseFaceIndex = %d, heBaseIndex = %d\n", baseFaceIndex, heBaseIndex);



	// assign new faces indices and new half edges indices
	std::vector<size_t> f;
	std::vector<size_t> nHeI; std::vector <half_edge> nHe;
	std::vector<size_t> nHeTwinI; std::vector <half_edge> nHeTwin;

	//printf("\nfaces: %d", m_faceData.size());
	for (int i = n - 3; i >= 0; i--)
	{
		f.push_back(baseFaceIndex + i);
	}

	//printf("n = %d", n);
	for (int i = 0; i < n - 3; i++)
	{
		nHeI.push_back(heBaseIndex + i);
		nHeTwinI.push_back(heBaseIndex + i + (n - 3));
	}

	//printf("f = ");
	//for (int i = 0; i < f.size(); i++)printf("%d ", f[i]);
	//printf("\nnHeI = ");
	//for (int i = 0; i < nHeI.size(); i++)printf("%d ", nHeI[i]);
	//printf("\nnHeTwinI = ");
	//for (int i = 0; i < nHeI.size(); i++)printf("%d ", nHeTwinI[i]);
	//printf("\n");



	// triangulating the hole

	m_heData[bIndex[0]].face = f[0];

	int i, j;
	i = j = -1;

	for (int x = 1; x < bIndex.size() - 2; x++)
	{
		i = n - x;
		j = x - 1;

		nHe.push_back(half_edge());
		nHeTwin.push_back(half_edge());

		m_heData[bIndex[i]].face = f[j];
		m_heData[bIndex[i]].next = nHeI[j];


		// half edge
		nHe[j].vert = v[j + 1];
		nHe[j].face = f[j];
		if (j == 0) nHe[j].next = bIndex[0];
		else nHe[j].next = nHeTwinI[j - 1];
		nHe[j].twin = nHeTwinI[j];

		//half edge twin
		nHeTwin[j].vert = ov;
		nHeTwin[j].face = f[j + 1];
		nHeTwin[j].next = bIndex[i - 1];
		nHeTwin[j].twin = nHeI[j];
	}

	if (i != -1 || j != -1){
		m_heData[bIndex[i - 1]].face = f[f.size() - 1];

		m_heData[bIndex[i - 2]].face = f[f.size() - 1];
		m_heData[bIndex[i - 2]].next = nHeTwinI[j];


		//printf("\n"); for (int i = 0; i < bIndex.size(); i++)printf("half edge [%d]: f:%d, n:%d, t:%d, v:%d\n", bIndex[i], m_heData[bIndex[i]].face, m_heData[bIndex[i]].next, m_heData[bIndex[i]].twin, m_heData[bIndex[i]].vert);
		//printf("\n"); for (int i = 0; i < nHe.size(); i++)printf("nHe [%d]: f:%d, n:%d, t:%d, v:%d\n", i, nHe[i].face, nHe[i].next, nHe[i].twin, nHe[i].vert);
		//printf("\n"); for (int i = 0; i < nHeTwin.size(); i++)printf("nHeTwin [%d]: f:%d, n:%d, t:%d, v:%d\n", i, nHeTwin[i].face, nHeTwin[i].next, nHeTwin[i].twin, nHeTwin[i].vert);

		for (int x = 0; x < nHe.size(); x++)
			m_heData.push_back(nHe[x]);

		for (int x = 0; x < nHeTwin.size(); x++)
			m_heData.push_back(nHeTwin[x]);
	}
	else
	{
		//printf("bIndex.size() = %d, f.size() = %d", bIndex.size(), f.size());
		m_heData[bIndex[1]].face = f[0];
		m_heData[bIndex[2]].face = f[0];
	}



	for (int x = 2; x < bIndex.size(); x++)
		m_faceData.push_back(bIndex[x]);

}

void EditMesh::distributeErrorToNeighbors(size_t vertIndex)
{
	std::vector<int> neighbors = getNeighborsIndices(vertIndex);
	float distributedError = 0.0f;

	for (int i = 0; i < vertices.size(); i++)
	if (vertices.at(i).vertexIndex == vertIndex)
	{
		distributedError = vertices.at(i).error / (float)neighbors.size();
		//set as HOLE_INDEX
		//vertices[i].vertexIndex = HOLE_INDEX;
		vertices[i].error = HOLE_INDEX;
		break;
	}

	for (int i = 0; i < neighbors.size(); i++)
	for (int j = 0; j < vertices.size(); j++)
	{
		if (vertices.at(j).vertexIndex == neighbors.at(i))
		{
			vertices.at(j).error += distributedError;
			break;
		}
	}

}

void EditMesh::sortVerticesByError()
{
	vertices_error v;

	int y = vertices.size() - 1;
	for (int x = 0; x<y; x++){
		for (int z = 0; z<y; z++){
			v.error = z;
			if (vertices[z].error < vertices[z + 1].error){
				v = vertices[z];
				vertices[z] = vertices[z + 1];
				vertices[z + 1] = v;
			}
		}
	}

}

float EditMesh::computeErrorAtVertex(size_t vertIndex){

	Eigen::Vector3d vect1;
	Eigen::Vector3d vect2;
	std::vector<float> alpha;
	float sumAlpha = 0;
	float error = 0;

	int vi = vertIndex;

	std::vector<Eigen::Vector3d> neighbors = getNeighborsXYZ(vi);

	//printf("\nvertex %d: (%1.2f, %1.2f,%1.2f)\n\n",vi, m_vertices[vi].x(), m_vertices[vi].y(), m_vertices[vi].z());

	for (int i = 0; i < neighbors.size(); i++)
	{
		//printf("neighbor %d: (%1.2f, %1.2f, %1.2f)\n", i, neighbors.at(i).x(), neighbors.at(i).y(), neighbors.at(i).z());

		if (i != neighbors.size() - 1){
			vect1 = neighbors.at(i) - m_vertices[vi];
			vect2 = neighbors.at(i + 1) - m_vertices[vi];
		}
		else
		{
			vect1 = neighbors.at(i) - m_vertices[vi];
			vect2 = neighbors.at(0) - m_vertices[vi];
		}

		alpha.push_back(acos((vect1.dot(vect2) / (sqrt(pow(vect1.x(), 2) + pow(vect1.y(), 2) + pow(vect1.z(), 2)) * sqrt(pow(vect2.x(), 2) + pow(vect2.y(), 2) + pow(vect2.z(), 2))))));

		//printf("alpha [%d] = %1.2f\n\n", i, alpha.at(i)*57.2957795);

		sumAlpha += alpha.at(i)*57.2957795;
	}

	error = sumAlpha / (2 * 3.14159265359);

	//printf("Error = %1.2f\n", error);

	return error;
}

std::vector<Eigen::Vector3d> EditMesh::getNeighborsXYZ(size_t vertIndex)
{
	std::vector<Eigen::Vector3d> neighbors;

	std::vector<half_edge> connectedHe;

	connectedHe = getConnectedHalfEdges(vertIndex);

	for (int i = 0; i < connectedHe.size(); i++)
	{
		neighbors.push_back(m_vertices[connectedHe.at(i).vert]);
	}

	return neighbors;
}

std::vector<int> EditMesh::getNeighborsIndices(size_t vertIndex)
{
	std::vector<int> neighbors;

	std::vector<half_edge> connectedHe;

	connectedHe = getConnectedHalfEdges(vertIndex);

	for (int i = 0; i < connectedHe.size(); i++)
	{
		neighbors.push_back(connectedHe.at(i).vert);
	}

	return neighbors;

}

std::vector<half_edge> EditMesh::getConnectedHalfEdges(size_t vertIndex)
{

	int firstHe, he;
	std::vector<half_edge> heVector;

	firstHe = m_heData[m_vertData[vertIndex]].twin;
	heVector.push_back(m_heData[firstHe]);

	he = firstHe;

	he = m_heData[m_heData[he].next].twin;

	while (he != firstHe){

		heVector.push_back(m_heData[he]);
		he = m_heData[m_heData[he].next].twin;
	}

	return heVector;
}

std::vector<size_t> EditMesh::getConnectedHalfEdgesIndices(size_t vertIndex)
{
	int firstHe, he;
	std::vector<size_t> heVector;

	firstHe = m_heData[m_vertData[vertIndex]].twin;
	heVector.push_back(firstHe);

	he = firstHe;

	he = m_heData[m_heData[he].next].twin;

	while (he != firstHe){

		heVector.push_back(he);
		he = m_heData[m_heData[he].next].twin;
	}

	return heVector;
}
#pragma endregion
// ------------------- Marching Cubes

void EditMesh::createImagesGrid(std::vector<Eigen::MatrixXd> images)
{
	scaleFactor = 0.01;
	//scaleFactor = 1;
	zScaleFactor = 0.75;
	//zScaleFactor = 1;

	size_t rows = images[1].innerSize();
	size_t cols = images[1].outerSize();
	size_t slices = images.size();

	std::cout << std::endl << std::endl;

	for (int i = 0; i < rows; i++)
		mainGrid.row.push_back(i);

	for (int j = 0; j < cols; j++)
		mainGrid.col.push_back(j);

	for (int k = 0; k < slices; k++)
		mainGrid.slice.push_back(k);

	//for (int i = 0; i < rows; i++)
	//{	
	//	for (int j = 0; j < cols; j++)
	//	{		
	//		for (int k = 0; k < images.size(); k++)
	//		{
	//			std::cout << "(" << mainGrid.row[i] << "," << mainGrid.col[j] << "," << mainGrid.slice[k] << ") , " << std::endl;
	//		}
	//	}
	//}

	mainGrid.isolevel = 0.5;

	std::cout << "row = " << mainGrid.row.size() << ", col = " << mainGrid.col.size() << ", slice = " << mainGrid.slice.size() << std::endl;

	/*mainGrid.gridCellMat = new grid_cell**[mainGrid.row.size()-1]();
	int row = mainGrid.row.size();
	int col = mainGrid.col.size();
	int slice = mainGrid.slice.size();
	mainGrid.gridCellMat = new grid_cell**[row-1]();
	for (int i = 0; i < row-1; ++i)
	{
		mainGrid.gridCellMat[i] = new grid_cell*[col-1]();
		for (int j = 0; j < col-1; j++)
		{
			mainGrid.gridCellMat[i][j] = new grid_cell[slice-1]();
			for (int k = 0; k < slice-1; k++)
				mainGrid.gridCellMat[i][j][k].val[0] = -1;
		}
	}*/

	initializeTables();
}

void EditMesh::createGridCells(std::vector<Eigen::MatrixXd> images)
{
	size_t HUGE_NUMBER = 3435973836;

	grid_cell iteratorCell;

	int cubeCount = 0;
	for (int k = 0; k < mainGrid.slice.size() - 1 ; k++)
	{
		std::cout << "\n============\nslice: " << k << std::endl;

		for (int j = 0; j < mainGrid.col.size() - 1; j++)
		{
			//std::cout << "row: " << j << std::endl;

			for (int i = 0; i < mainGrid.row.size() - 1; i++){

				//std::cout << "col: " << i << std::endl;

				//std::cout << "cube # " << cubeCount << std::endl;
				iteratorCell.cellIndex = cubeCount;

#pragma region creating grid cells

				iteratorCell.pts[0] = Eigen::Vector3d(mainGrid.row[0 + i], mainGrid.col[0 + j], mainGrid.slice[0 + k]);
				iteratorCell.pts[1] = Eigen::Vector3d(mainGrid.row[0 + i], mainGrid.col[1 + j], mainGrid.slice[0 + k]);
				iteratorCell.pts[2] = Eigen::Vector3d(mainGrid.row[1 + i], mainGrid.col[1 + j], mainGrid.slice[0 + k]);
				iteratorCell.pts[3] = Eigen::Vector3d(mainGrid.row[1 + i], mainGrid.col[0 + j], mainGrid.slice[0 + k]);

				iteratorCell.pts[4] = Eigen::Vector3d(mainGrid.row[0 + i], mainGrid.col[0 + j], mainGrid.slice[1 + k]);
				iteratorCell.pts[5] = Eigen::Vector3d(mainGrid.row[0 + i], mainGrid.col[1 + j], mainGrid.slice[1 + k]);
				iteratorCell.pts[6] = Eigen::Vector3d(mainGrid.row[1 + i], mainGrid.col[1 + j], mainGrid.slice[1 + k]);
				iteratorCell.pts[7] = Eigen::Vector3d(mainGrid.row[1 + i], mainGrid.col[0 + j], mainGrid.slice[1 + k]);

				iteratorCell.val[0] = images[0 + k](0 + i, 0 + j);
				iteratorCell.val[1] = images[0 + k](0 + i, 1 + j);
				iteratorCell.val[2] = images[0 + k](1 + i, 1 + j);
				iteratorCell.val[3] = images[0 + k](1 + i, 0 + j);

				iteratorCell.val[4] = images[1 + k](0 + i, 0 + j);
				iteratorCell.val[5] = images[1 + k](0 + i, 1 + j);
				iteratorCell.val[6] = images[1 + k](1 + i, 1 + j);
				iteratorCell.val[7] = images[1 + k](1 + i, 0 + j);


				iteratorCell.vertices[0] = (iteratorCell.pts[0] + iteratorCell.pts[1]) / 2;
				iteratorCell.vertices[1] = (iteratorCell.pts[1] + iteratorCell.pts[2]) / 2;
				iteratorCell.vertices[2] = (iteratorCell.pts[2] + iteratorCell.pts[3]) / 2;
				iteratorCell.vertices[3] = (iteratorCell.pts[3] + iteratorCell.pts[0]) / 2;

				iteratorCell.vertices[4] = (iteratorCell.pts[4] + iteratorCell.pts[5]) / 2;
				iteratorCell.vertices[5] = (iteratorCell.pts[5] + iteratorCell.pts[6]) / 2;
				iteratorCell.vertices[6] = (iteratorCell.pts[6] + iteratorCell.pts[7]) / 2;
				iteratorCell.vertices[7] = (iteratorCell.pts[7] + iteratorCell.pts[4]) / 2;

				iteratorCell.vertices[8] = (iteratorCell.pts[4] + iteratorCell.pts[0]) / 2;
				iteratorCell.vertices[9] = (iteratorCell.pts[5] + iteratorCell.pts[1]) / 2;
				iteratorCell.vertices[10] = (iteratorCell.pts[6] + iteratorCell.pts[2]) / 2;
				iteratorCell.vertices[11] = (iteratorCell.pts[7] + iteratorCell.pts[3]) / 2;

				for (int h = 0; h < 12; h++)
					iteratorCell.vertices_global[h] = HUGE_NUMBER;

				//for (int i = 0; i < 8; i++)
				//{
				//	std::cout << "pt[" << i << "]: (" << iteratorCell.pts[i].x()<< "," << iteratorCell.pts[i].y() <<","<< iteratorCell.pts[i].z() << ") = ";
				//	std::cout << iteratorCell.val[i] << std::endl;
				//}

				//for (int i = 0; i < 12; i++)
				//{
				//	std::cout << "vertices[" << i << "]: (" << iteratorCell.vertices[i].x() << "," << iteratorCell.vertices[i].y() << "," << iteratorCell.vertices[i].z() << ") "<<std::endl;
				//}


#pragma endregion

#pragma region computing cell index

				iteratorCell.cubeindex = 0;
				if (iteratorCell.val[0] < mainGrid.isolevel) iteratorCell.cubeindex |= 1;
				if (iteratorCell.val[1] < mainGrid.isolevel) iteratorCell.cubeindex |= 2;
				if (iteratorCell.val[2] < mainGrid.isolevel) iteratorCell.cubeindex |= 4;
				if (iteratorCell.val[3] < mainGrid.isolevel) iteratorCell.cubeindex |= 8;
				if (iteratorCell.val[4] < mainGrid.isolevel) iteratorCell.cubeindex |= 16;
				if (iteratorCell.val[5] < mainGrid.isolevel) iteratorCell.cubeindex |= 32;
				if (iteratorCell.val[6] < mainGrid.isolevel) iteratorCell.cubeindex |= 64;
				if (iteratorCell.val[7] < mainGrid.isolevel) iteratorCell.cubeindex |= 128;

				//std::cout << "cubeIndex = ";
				//printf("%d\n", iteratorCell.cubeindex);

				//iteratorCell.edgeTableVal = edgeTable[iteratorCell.cubeindex];
				//std::cout << "EdgeTableVal = " << iteratorCell.edgeTableVal << std::endl;

				//std::cout << "EdgeTableVal = {";

				Eigen::Vector3d currentVert;
				size_t vertIndex;
				for (int c = 0; c < 16; c++)
				{
					iteratorCell.triTableVal[c] = triTable[iteratorCell.cubeindex][c];
					//std::cout << iteratorCell.triTableVal[c] << ", ";

					if (iteratorCell.triTableVal[c] != -1)
					{

						currentVert = iteratorCell.vertices[iteratorCell.triTableVal[c]];
						vertIndex = getVertexIndex(currentVert);
						//vertIndex = getGlobalVertexIndex(currentVert, cubeCount);
						iteratorCell.vertices_global[iteratorCell.triTableVal[c]] = vertIndex;

						triangles.push_back(vertIndex);
						//std::cout << "[" << vertIndex << "]: (" << currentVert.x() << "," << currentVert.y() << "," << currentVert.z() << ")" << std::endl;
					}
				}

				

				/*for (int h = 0; h < 12; h++)
				{
					std::cout << "vertGlobal[" << h << "]: " << iteratorCell.vertices_global[h] << std::endl;
				}*/

				//std::cout << "}" << std::endl;

				//std::cout << "(" << i << ", " << j << ", " << k << "): #[" << iteratorCell.cellIndex << "]: cubeIndex = " << iteratorCell.cubeindex << std::endl;


#pragma endregion

				mainGrid.gridCells.push_back(iteratorCell);
				//mainGrid.gridCellMat[i][j][k] = iteratorCell;
				cubeCount++;
			}
		}
	}

	for (int i = 0; i < verts.size(); i++)
	{
		vertsInput.push_back(verts[i].x()*scaleFactor);
		vertsInput.push_back(verts[i].y()*scaleFactor);
		vertsInput.push_back(verts[i].z()*scaleFactor*zScaleFactor);
	}

	//std::cout << "Vertices: " << std::endl;
	//for (int i = 0; i < vertsInput.size(); i++)
	//{
	//	std::cout << vertsInput[i] << std::endl;
	//}

	//std::cout << "\n===========\nTriangles:" << std::endl;
	//for (int i = 0; i < triangles.size(); i++)
	//{
	//	std::cout << triangles[i] << std::endl;
	//}

	this->init(vertsInput, triangles);
	edit_count++;

}

size_t EditMesh::getGlobalVertexIndex(Eigen::Vector3d vertex, size_t cubeCount)
{
	size_t HUGE_NUMBER = 3435973836;
	size_t index = -1;
	size_t count = 0;
	grid_cell cell;

	for (int i = 0; i < mainGrid.slice.size() - 1; i++)
	for (int j = 0; j < mainGrid.col.size() - 1; j++)
	for (int k = 0; k < mainGrid.row.size() - 1; k++)
	{
		count++;
		if (count >= cubeCount)
			break;
		else
		{
			cell = mainGrid.gridCellMat[i][j][k];
			for (int c = 0; c < 12; c++)
			if (cell.vertices_global[c]!=HUGE_NUMBER)
			if (vertex.x() == cell.vertices[c].x() && vertex.y() == cell.vertices[c].y() && vertex.z() == cell.vertices[c].z())
			{
				index = cell.vertices_global[c];
				break;
			}
		}
	}


	if (index == -1)
	{
		verts.push_back(vertex);
		index = verts.size() - 1;
	}

	return index;
}

size_t EditMesh::getVertexIndex(Eigen::Vector3d vertex)
{
	size_t index = -1;
	for (int i = 0; i < verts.size(); i++)
	{
		if (vertex.x() == verts[i].x() && vertex.y() == verts[i].y() && vertex.z() == verts[i].z())
		{
			index = i;
			break;
		}
	}

	if (index == -1)
	{
		verts.push_back(vertex);
		index = verts.size() - 1;
	}

	return index;
}

void EditMesh::initializeTables()
{
	int edgeTableLocal[256] = {
		0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
		0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
		0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
		0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
		0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
		0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
		0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
		0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
		0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
		0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
		0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
		0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
		0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
		0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
		0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
		0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
		0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
		0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
		0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
		0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
		0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
		0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
		0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
		0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
		0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
		0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
		0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
		0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
		0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
		0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0 };

	for (int i = 0; i < 256; i++)
	{
		edgeTable[i] = edgeTableLocal[i];
	}

	int triTableLocal[256][16] =
	{ { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
	{ 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
	{ 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
	{ 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
	{ 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
	{ 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
	{ 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
	{ 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
	{ 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
	{ 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
	{ 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
	{ 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
	{ 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
	{ 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
	{ 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
	{ 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
	{ 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
	{ 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
	{ 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
	{ 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
	{ 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
	{ 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
	{ 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
	{ 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
	{ 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
	{ 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
	{ 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
	{ 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
	{ 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
	{ 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
	{ 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
	{ 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
	{ 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
	{ 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
	{ 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
	{ 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
	{ 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
	{ 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
	{ 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
	{ 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
	{ 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
	{ 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
	{ 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
	{ 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
	{ 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
	{ 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
	{ 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
	{ 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1 },
	{ 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1 },
	{ 10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1 },
	{ 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1 },
	{ 7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1 },
	{ 6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1 },
	{ 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1 },
	{ 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1 },
	{ 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1 },
	{ 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1 },
	{ 10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1 },
	{ 10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1 },
	{ 9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1 },
	{ 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1 },
	{ 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1 },
	{ 7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1 },
	{ 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1 },
	{ 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1 },
	{ 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1 },
	{ 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1 },
	{ 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1 },
	{ 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1 },
	{ 6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1 },
	{ 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1 },
	{ 6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1 },
	{ 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1 },
	{ 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1 },
	{ 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1 },
	{ 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1 },
	{ 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1 },
	{ 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1 },
	{ 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1 },
	{ 0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1 },
	{ 11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1 },
	{ 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1 },
	{ 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1 },
	{ 2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1 },
	{ 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1 },
	{ 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1 },
	{ 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1 },
	{ 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1 },
	{ 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1 },
	{ 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1 },
	{ 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1 },
	{ 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1 },
	{ 9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1 },
	{ 5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1 },
	{ 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1 },
	{ 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1 },
	{ 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1 },
	{ 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1 },
	{ 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1 },
	{ 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1 },
	{ 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1 },
	{ 11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1 },
	{ 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1 },
	{ 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1 },
	{ 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1 },
	{ 1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1 },
	{ 4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1 },
	{ 0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1 },
	{ 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } };

	for (int i = 0; i < 256; i++)
	for (int j = 0; j < 16; j++)
		triTable[i][j] = triTableLocal[i][j];
}