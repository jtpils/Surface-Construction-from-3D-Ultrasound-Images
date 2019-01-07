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
#ifndef OBJECT_MESH_H
#define OBJECT_MESH_H

#include <iostream>
#include "DrawMesh.h"
#include "EditMesh.h"
#include <memory>

typedef std::shared_ptr<EditMesh> EditMesh_ptr;
typedef std::shared_ptr<DrawMesh> DrawMesh_ptr;

class Model
{
public:
    Model();
    Model::Model(EditMesh_ptr em);
    Model::~Model();

    void init(RenderState &state);

    void drawMesh(int primitive = GL_TRIANGLES);

    /***********************************************************
     * the reflection interface to get mesh information
     ***********************************************************/
    std::size_t  info_sizev() {return m_em->get_vert_size();}
    std::size_t  info_sizef() {return m_em->get_face_size();}
	Eigen::Vector3d info_vertex(std::size_t i);
	Eigen::Vector3d info_vertexNormal(std::size_t i);

    void info_bbox(Eigen::Vector3d &bboxMin, Eigen::Vector3d &bboxMax);
	void getIndicesForFace(size_t tri_index, size_t indicesForFace[3]);

    /***********************************************
     * mesh modification algorithms
     ***********************************************/
    // expose your editmesh functions here
	bool edit_subdivision(std::size_t subdivMethod = 0);
	bool edit_set_num_simplify(std::size_t numVertRemove = 1);
	bool edit_set_error_simplify(std::size_t eorrorRemove = 1);
	void edit_forward_simplify();
	void edit_backward_simplify();
	/***********************************************
     * re_meshing algorithms
     ***********************************************/
	void edit_remesh();
	void edit_loadMesh();
	void edit_load30Slices();
	void edit_load56Slices();
	void edit_load113Slices();
	///////////////////////////////////////////////////////////////////////////////////
	const EditMesh_ptr get_editMesh() const;
    void               set_editMesh(EditMesh_ptr em);

protected:
    void updateDrawMesh();

    /***********************************************
     * Variables for this Instance of the EditMesh
     ***********************************************/
    // the bounding box for this instance of the geometry
    Eigen::Vector3f bboxMin;
    Eigen::Vector3f bboxMax;

    //indicates the last drawn edit mesh, if this does not match
    // the value in m_em then the draw mesh needs to be updated
    int last_drawn;

    DrawMesh_ptr m_dm;
    EditMesh_ptr m_em;
};

inline bool Model::edit_subdivision(std::size_t subdivMethod)
{
	if (m_em == NULL) 
		{
			std::cerr << "there is no mesh to be subdivided!" << std::endl;
			return false;
		}

	m_em->Subdivision(subdivMethod); // subdivMethod: 0 for butterfly, 1 for loop;
	return true;
}

inline bool Model::edit_set_num_simplify(std::size_t numVertRemove)
{
	if (m_em == NULL) 
	{
		std::cerr << "there is no mesh to be simplified!" << std::endl;
		return false;
	}

	m_em->num_vert_remove(numVertRemove); // set the number of vertices to be removed
	return true;
}

inline bool Model::edit_set_error_simplify(std::size_t err)
{
	if (m_em == NULL) 
	{
		std::cerr << "there is no mesh to be simplified!" << std::endl;
		return false;
	}

	m_em->Errormetric_simplify(err); // set the number of vertices to be removed
	return true;
}

inline void Model::edit_forward_simplify()
{
	m_em->oneStep_simplify();
}

inline void Model::edit_backward_simplify()
{
	m_em->mesh_backward_simplify(); 
}
inline void Model::edit_remesh(){
	//m_em->refine_mesh();
	m_em->re_mesh();
}
inline void Model::edit_loadMesh(){
	m_em->loadMesh();
}
inline void Model::edit_load30Slices(){m_em->load30Slices();}
inline void Model::edit_load56Slices(){m_em->load56Slices();}
inline void Model::edit_load113Slices(){m_em->load113Slices();}

inline const EditMesh_ptr Model::get_editMesh() const {
	return m_em;
}

inline void Model::set_editMesh(EditMesh_ptr em) {
    m_em = em;
}

#endif //OBJECT_MESH_H