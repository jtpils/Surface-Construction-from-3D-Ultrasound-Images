/* Copyright (c) Russell Gillette
 * December 2013
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

/* == MeshUtils.h ==
 * 
 * Utility and helper functions that operate on meshes or mesh data without a
 * specified Mesh instance.
 */

#ifndef MESH_UTILS_H
#define MESH_UTILS_H

#include "Model.h"
#include "DrawMesh.h"
#include "EditMesh.h"
#include "OBJLoader.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Eigen/Dense"

#include <string>
 
#define DEFAULT_END_FRAME 1000000
const double PI=3.1415926535897;


// load an Edit Mesh object from an OBJ file.
EditMesh *loadEditMeshFromFile(string file_name)
{
    EditMesh *m = NULL;

    // parse mesh
    ObjLoader obj_parser(file_name);
    //obj_parser.generateNormals(); // don't use normals atm so no point in doing this

    // gather data from the parser to construct mesh
    GLubyte *data = NULL;
    int *indices;
    int data_size, num_indices, num_attr;
    obj_attrib_info *parsed_attr_info;

    // export parser data
    // NOTE: data is owned by the obj_parser
    obj_parser.objExportGLSeparate(data_size, data, num_indices, indices, num_attr, parsed_attr_info);

    if (parsed_attr_info == NULL)
    {
        printf("Mesh Not Found: Failed to load\n");
        return NULL;
    }

    std::vector<double> xyzPositions;
    std::vector<std::size_t> triangleIndices;

    // populate indices
    for (int i = 0; i < num_indices; i++)
        triangleIndices.push_back(indices[i]);
    
    // populate vertices from byte array (so lots of pointer pushing)
    int attr_end = data_size;
    int v_stride = parsed_attr_info[0].data_stride;
    int v_offset = parsed_attr_info[0].data_offset;

    if (v_stride == 0)
    for (int i = 0; i < num_attr; i++)
    {
        int off = parsed_attr_info[i].data_offset;
        if (off < attr_end && off > v_offset)
            attr_end = off;
    }

    int      attrib_size = parsed_attr_info[0].attrib_size;
    int      elem_size = attrib_size / parsed_attr_info[0].num_comp;

    GLubyte *pAttrib = data;
    GLubyte *pEnd    = data + attr_end;

    // TODO: safety check on number of elements per attribute
    // (there should never be less than 3 but who knows)
    // if (sizeof(float) != elem_size)
    //    unhandled as of right now

    // not particularily safe...
    for (; pAttrib < pEnd; pAttrib += elem_size)
    {
        double tmp = (double)*(float*)pAttrib;
        xyzPositions.push_back(tmp);
    }

    m = new EditMesh();
    m->init(xyzPositions, triangleIndices);

    return m;
}

// load a Mesh object from an OBJ file.
// NOTE: this mesh is packaged for drawing, and thus difficult to modify
//    Mesh for editable mesh
DrawMesh *loadDrawMeshFromFile(RenderState &state, string file_name)
{
    DrawMesh *m = NULL;

    // parse mesh
    ObjLoader obj_parser(file_name);
    obj_parser.generateNormals();

    // gather data from the parser to construct mesh
    GLubyte *data = NULL;
    int *indices;
    int data_size, num_indices, num_attr;
    obj_attrib_info *parsed_attr_info;
    attrib_info     *attr_info;

    // export parser data
    // NOTE: data is owned by the obj_parser
    obj_parser.objExportGLSeparate(data_size, data, num_indices, indices, num_attr, parsed_attr_info);

    attr_info = new attrib_info[num_attr];
    for (int j = 0; j < num_attr; j++)
    {
        attr_info[j].attrib_number = parsed_attr_info[j].attrib_number;
        attr_info[j].attrib_size   = parsed_attr_info[j].attrib_size;
        attr_info[j].data_offset   = parsed_attr_info[j].data_offset;
        attr_info[j].data_stride   = parsed_attr_info[j].data_stride;
        attr_info[j].num_comp      = parsed_attr_info[j].num_comp;
    }

    m = new DrawMesh(state);
    m->init(1);
    m->loadVBuffer(0, data_size, data, 0, num_attr, attr_info);
    m->loadIBuffer(num_indices, sizeof(int), indices);

    delete [] attr_info;

    return m;
}

Model *loadModelFromFile(RenderState &state, string file_name)
{
	EditMesh *em = loadEditMeshFromFile(file_name);
    Model *m = new Model(EditMesh_ptr(em));
    m->init(state);
    return m;
}

DrawMesh *createAxis(RenderState & state, float scale)
{
    DrawMesh *m;
    attrib_info details[2];

    // NOTE: vertices are duplicated, as theu have different normals
    // and this is not possible to represent in opengl
    float data[] = 
    {
        // vertices
         0.0f,  1.0f,  0.0f,  -0.04f, 0.0f,  0.04f,  0.04f, 0.0f,  0.04f,
         0.0f,  1.0f,  0.0f,  -0.04f, 0.0f, -0.04f, -0.04f, 0.0f,  0.04f,
         0.0f,  1.0f,  0.0f,   0.04f, 0.0f, -0.04f, -0.04f, 0.0f, -0.04f,
         0.0f,  1.0f,  0.0f,   0.04f, 0.0f,  0.04f,  0.04f, 0.0f, -0.04f,
         0.04f, 0.0f,  0.04f,  0.04f, 0.0f, -0.04f, -0.04f, 0.0f, -0.04f,
        -0.04f, 0.0f, -0.04f, -0.04f, 0.0f,  0.04f,  0.04f, 0.0f,  0.04f,
        // start normals: this is a hack to approximate reasonable normals
         0.0f,  0.0f,  1.0f,   0.0f, 0.0f,   1.0f,   0.0f,  0.0f,  1.0f,
        -1.0f,  0.0f,  0.0f,  -1.0f, 0.0f,   0.0f,  -1.0f,  0.0f,  0.0f,
         0.0f,  0.0f, -1.0f,   0.0f, 0.0f,  -1.0f,   0.0f,  0.0f, -1.0f,
         1.0f,  0.0f,  0.0f,   1.0f, 0.0f,   0.0f,   1.0f,  0.0f,  0.0f,
         0.0f, -1.0f,  0.0f,   0.0f, -1.0f,  0.0f,   0.0f, -1.0f,  0.0f,
         0.0f, -1.0f,  0.0f,   0.0f, -1.0f,  0.0f,   0.0f, -1.0f,  0.0f 
    };

    // scale the vertex locations relative to the scale
    // parameter (note I only go to 54, the num verts)
    for (int i = 0; i < 54; i++)
    {
        data[i] *= scale;
    }

    // since we have redefined our vertices, the indices
    // are just a linear ordering
    int indices[18] = {0,1,2,3,4,5,6,7,8,9,10,11,12,
        13,14,15,16,17};

    // attribute information for vertices
    details[0].attrib_number = 0;
    details[0].attrib_size   = 3 * sizeof(float); // size of one vertex
    details[0].data_offset   = 0; // byte offset to vertices
    details[0].data_stride   = 0; // space between vertices
    details[0].num_comp      = 3;

    // attribute information for normals
    details[1].attrib_number = 1;
    details[1].attrib_size   = 3 * sizeof(float); // size of one normal
    details[1].data_offset   = 18 * 3 * sizeof(float); // byte offset to normals
    details[1].data_stride   = 0; // space between normals
    details[1].num_comp      = 3;

    m = new DrawMesh(state);
    m->init(1);
    m->loadVBuffer(0, sizeof(data), (GLubyte*)data, 0, 2, details);
    m->loadIBuffer(18, sizeof(int), indices);

    return m;
}

DrawMesh *createGem(RenderState & state, float scale)
{
    DrawMesh *m;
    attrib_info details[2];

    // NOTE: vertices are duplicated, as they have different normals
    // and this is not possible to represent in opengl
    float data[] = 
    {
         // vertices
         0.5f, -0.5f,  0.0f, //v1
         0.5f,  0.5f,  0.0f, //v2
         0.0f,  0.0f,  1.0f, //v5

         0.5f,  0.5f,  0.0f, //v2
        -0.5f,  0.5f,  0.0f, //v3
         0.0f,  0.0f,  1.0f, //v5

        -0.5f,  0.5f,  0.0f, //v3
        -0.5f, -0.5f,  0.0f, //v4
         0.0f,  0.0f,  1.0f, //v5

        -0.5f, -0.5f,  0.0f, //v4
         0.5f, -0.5f,  0.0f, //v1
         0.0f,  0.0f,  1.0f, //v5

         0.5f,  0.5f,  0.0f, //v2
         0.5f, -0.5f,  0.0f, //v1
         0.0f,  0.0f, -1.0f, //v6

        -0.5f,  0.5f,  0.0f, //v3
         0.5f,  0.5f,  0.0f, //v2
         0.0f,  0.0f, -1.0f, //v6

        -0.5f, -0.5f,  0.0f, //v4
        -0.5f,  0.5f,  0.0f, //v3
         0.0f,  0.0f, -1.0f, //v6

         0.5f, -0.5f,  0.0f, //v1
        -0.5f, -0.5f,  0.0f, //v4
         0.0f,  0.0f, -1.0f, //v6

         // normals
         0.894427191f,  0.0f,  0.447213595f, //n1
         0.894427191f,  0.0f,  0.447213595f, //n1
         0.894427191f,  0.0f,  0.447213595f, //n1

         0.0f,  0.894427191f,  0.447213595f, //n2
         0.0f,  0.894427191f,  0.447213595f, //n2
         0.0f,  0.894427191f,  0.447213595f, //n2

        -0.894427191f,  0.0f,  0.447213595f, //n3
        -0.894427191f,  0.0f,  0.447213595f, //n3
        -0.894427191f,  0.0f,  0.447213595f, //n3

         0.0f, -0.894427191f,  0.447213595f, //n4
         0.0f, -0.894427191f,  0.447213595f, //n4
         0.0f, -0.894427191f,  0.447213595f, //n4

         0.894427191f,  0.0f, -0.447213595f, //n5
         0.894427191f,  0.0f, -0.447213595f, //n5
         0.894427191f,  0.0f, -0.447213595f, //n5

         0.0f,  0.894427191f, -0.447213595f, //n6
         0.0f,  0.894427191f, -0.447213595f, //n6
         0.0f,  0.894427191f, -0.447213595f, //n6

        -0.894427191f,  0.0f, -0.447213595f, //n7
        -0.894427191f,  0.0f, -0.447213595f, //n7
        -0.894427191f,  0.0f, -0.447213595f, //n7

         0.0f, -0.894427191f, -0.447213595f, //n8
         0.0f, -0.894427191f, -0.447213595f, //n8
         0.0f, -0.894427191f, -0.447213595f, //n8
    };

    // scale the vertex locations relative to the scale
    // parameter (note I only go to 72, the num verts)
    for (int i = 0; i < 72; i++)
    {
        data[i] *= scale;
    }

    // since we have redefined our vertices, the indices
    // are just a linear ordering
    int indices[24] = {0,1,2,3,4,5,6,7,8,9,10,11,12,
        13,14,15,16,17,18,19,20,21,22,23};

    // attribute information for vertices
    details[0].attrib_number = 0;
    details[0].attrib_size   = 3 * sizeof(float); // size of one vertex
    details[0].data_offset   = 0; // byte offset to vertices
    details[0].data_stride   = 0; // space between vertices
    details[0].num_comp      = 3;

    // attribute information for normals
    details[1].attrib_number = 1;
    details[1].attrib_size   = 3 * sizeof(float); // size of one normal
    details[1].data_offset   = 24 * 3 * sizeof(float); // byte offset to normals
    details[1].data_stride   = 0; // space between normals
    details[1].num_comp      = 3;

    m = new DrawMesh(state);
    m->init(1);
    m->loadVBuffer(0, sizeof(data), (GLubyte*)data, 0, 2, details);
    m->loadIBuffer(24, sizeof(int), indices);

    return m;
}

void drawBox(int x1, int y1, int x2, int y2)
{
    // box has zero area
    if (x1 == x2 && y1 == y2)
        return;

    // safety check to avoid triangle flipping
    if (x1 > x2)
    {
        int tmp = x2;
        x2 = x1;
        x1 = tmp;
    }
    if (y1 > y2)
    {
        int tmp = y2;
        y2 = y1;
        y1 = tmp;
    }
    float fx1 = ((float)x1/c_state.width)*2 - 1;
    float fx2 = ((float)x2/c_state.width)*2 - 1;
    float fy1 = ((float)y1/c_state.height)*2 - 1;
    float fy2 = ((float)y2/c_state.height)*2 - 1;

    GLfloat vertices[] = { fx1, fy1, 0.5f,
                           fx2, fy1, 0.5f,
                           fx1, fy2, 0.5f,
                           fx1, fy2, 0.5f,
                           fx2, fy1, 0.5f,
                           fx2, fy2, 0.5f };

    // activate and specify pointer to vertex array
    // this could be much more efficient if I wasn't lazy
    GLuint tmp;
    glGenBuffers(1, &tmp);
    glBindBuffer(GL_ARRAY_BUFFER, tmp);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*18, vertices, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    // draw box
    glDrawArrays(GL_TRIANGLES, 0, 18);

    // deactivate vertex arrays after drawing
    glDisableVertexAttribArray(0);
    glDeleteBuffers(1, &tmp);
}

//points is {vertex, endpoint, vertex, endpoint, ...} where each sequential
//vertex and endpoint define a line to be drawn
//colours is {colour, colour, colour, ...}, the colour of each point
//assumed same size as points
void drawLines(std::vector<glm::vec3> points, std::vector<glm::vec3> colour_vecs)
{
	int size = points.size();
	GLfloat *vertices = new GLfloat[3*size];
	GLfloat *colours  = new GLfloat[3*size];
	for (int i = 0; i < size; i++)
	{
		vertices[3*i]   = points[i].x;
		vertices[3*i+1] = points[i].y;
		vertices[3*i+2] = points[i].z;
		colours[3*i]    = colour_vecs[i].r;
		colours[3*i+1]  = colour_vecs[i].g;
		colours[3*i+2]  = colour_vecs[i].b;
	}

	// activate and specify pointer to vertex array
    // this could be much more efficient if I wasn't lazy
    GLuint tmp;
    glGenBuffers(1, &tmp);
    glBindBuffer(GL_ARRAY_BUFFER, tmp);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*size, vertices, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	GLuint tmp2;
	glGenBuffers(1, &tmp2);
	glBindBuffer(GL_ARRAY_BUFFER, tmp2);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*size, colours, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//attrib_number, num_comp, GL_FLOAT, GL_FALSE, data_stride, (GLvoid *) data_offset

    // draw the lines
	glLineWidth(2.0);
	//glEnable(GL_LINE_SMOOTH);
	glDrawArrays(GL_LINES, 0, 3*size);
	//glDisable(GL_LINE_SMOOTH);
	//glDrawArrays(GL_TRIANGLES, 0,  3*size);

    // deactivate vertex arrays after drawing
    glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
    glDeleteBuffers(1, &tmp);
	glDeleteBuffers(1, &tmp2);
	delete[] vertices;
	delete[] colours;

}

// returns a color indicating which quadtrant the modelview tranform has
// translated the point (0,0,0) to.
glm::vec3 identifyQuadrant(glm::mat4 modelview)
{
    glm::vec3 out = glm::vec3();
    out.r = (modelview[3].r > 0) ? 1 : 0;
    out.g = (modelview[3].g > 0) ? 1 : 0;
    out.b = (modelview[3].b > 0) ? 1 : 0;
    out *= 0.7;
    out += glm::vec3(0.3f);
    return out;
}



DrawMesh *createPainter(RenderState & state, float scale, EditMesh& editMesh, std::vector<size_t>& vertVecColor)
{
	DrawMesh *m;
	attrib_info details[2];

	// extract the vertices' position and face normal data
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertPositon;
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> faceNormal;
	for (size_t i(0), iEnd(vertVecColor.size()); i < iEnd; ++i)
	{
		vface_iterator iter; // initial the iterator visiting the one ring of "vertex"
		if (!editMesh.init_iterator(iter, vertVecColor[i]))
			std::cerr << "error because of float point" << std::endl;
		do {
			faceNormal.push_back(editMesh.get_normal(iter));

			vertPositon.push_back(editMesh.get_vertex(iter.m_cur->vert));
			vertPositon.push_back(editMesh.get_vertex(editMesh.next(*iter.m_cur).vert));
			vertPositon.push_back(editMesh.get_vertex(editMesh.prev(*iter.m_cur).vert));
	
		} while(editMesh.advance_iterator(iter));
	}

	assert (faceNormal.size() == vertPositon.size() / 3 && "number of face and vertex don't match");

	// Copy from faceNormal and vertPostion to data
	size_t num(faceNormal.size());
	float* data = new float[faceNormal.size() * 18];
	for (size_t i(0), iEnd(faceNormal.size()); i < iEnd; ++i)
	{
		data[i * 9 + 0] = vertPositon[i * 3 + 0][0];
		data[i * 9 + 1] = vertPositon[i * 3 + 0][1];
		data[i * 9 + 2] = vertPositon[i * 3 + 0][2];
		data[i * 9 + 3] = vertPositon[i * 3 + 1][0];
		data[i * 9 + 4] = vertPositon[i * 3 + 1][1];
		data[i * 9 + 5] = vertPositon[i * 3 + 1][2];
		data[i * 9 + 6] = vertPositon[i * 3 + 2][0];
		data[i * 9 + 7] = vertPositon[i * 3 + 2][1];
		data[i * 9 + 8] = vertPositon[i * 3 + 2][2];

		data[9 * faceNormal.size() + i * 9 + 0] = faceNormal[i][0];
		data[9 * faceNormal.size() + i * 9 + 1] = faceNormal[i][1];
		data[9 * faceNormal.size() + i * 9 + 2] = faceNormal[i][2];
		data[9 * faceNormal.size() + i * 9 + 3] = faceNormal[i][0];
		data[9 * faceNormal.size() + i * 9 + 4] = faceNormal[i][1];
		data[9 * faceNormal.size() + i * 9 + 5] = faceNormal[i][2];
		data[9 * faceNormal.size() + i * 9 + 6] = faceNormal[i][0];
		data[9 * faceNormal.size() + i * 9 + 7] = faceNormal[i][1];
		data[9 * faceNormal.size() + i * 9 + 8] = faceNormal[i][2];
	}


	// scale the vertex locations relative to the scale
	for (int i = 0; i < faceNormal.size() * 9; i++)
	{
		data[i] *= scale;
	}

	// since we have redefined our vertices, the indices
	// are just a linear ordering
	int* indices = new int[faceNormal.size() * 3];
	for (size_t i(0), iEnd(faceNormal.size() * 3); i < iEnd; ++i)
	{
		indices[i] = i;
	}
		
	// attribute information for vertices
	details[0].attrib_number = 0;
	details[0].attrib_size   = 3 * sizeof(float); // size of one vertex
	details[0].data_offset   = 0; // byte offset to vertices
	details[0].data_stride   = 0; // space between vertices
	details[0].num_comp      = 3;

	// attribute information for normals
	details[1].attrib_number = 1;
	details[1].attrib_size   = 3 * sizeof(float); // size of one normal
	details[1].data_offset   = faceNormal.size() * 3 * 3 * sizeof(float); // byte offset to normals
	details[1].data_stride   = 0; // space between normals
	details[1].num_comp      = 3;

	m = new DrawMesh(state);
	m->init(1);
	m->loadVBuffer(0, sizeof(data), (GLubyte*)data, 0, 2, details);
	m->loadIBuffer(faceNormal.size() * 3, sizeof(int), indices);

	delete[] data;
	delete[] indices;

	return m;
}

void drawTriangles(EditMesh& editMesh, std::vector<size_t>& vertVecColor)
{	
	// extract the vertices' position and face normal data
	std::vector<glm::vec3> point_vecs, colour_vecs;
	std::vector<bool> isColored(editMesh.get_face_size(), false);
	for (size_t i(0), iEnd(vertVecColor.size()); i < iEnd; ++i)
	{
		vface_iterator iter; // initial the iterator visiting the one ring of "vertex"
		if (!editMesh.init_iterator(iter, vertVecColor[i]))
			std::cerr << "error because of float point" << std::endl;
		do {
			//if (i == 0)
			//	{
			//		colour_vecs.push_back(glm::vec3(1, 0, 0));
			//		colour_vecs.push_back(glm::vec3(1, 0, 0));
			//		//colour_vecs.push_back(glm::vec3(1, 0, 0));
			//	}
			//else
			//	{
			//		colour_vecs.push_back(glm::vec3(0, 1, 0));
			//		colour_vecs.push_back(glm::vec3(0, 1, 0));
			//		//colour_vecs.push_back(glm::vec3(0, 1, 0));
			//	}

			if (isColored[iter.m_cur->face] == true)
				continue;		
			
			Eigen::Vector3d v1(editMesh.get_vertex(iter.m_cur->vert));
			Eigen::Vector3d v2(editMesh.get_vertex(editMesh.next(*iter.m_cur).vert));
			Eigen::Vector3d v3(editMesh.get_vertex(editMesh.prev(*iter.m_cur).vert));

			// offset the vertices a bit
			Eigen::Vector3d n1(editMesh.get_vnormal(iter.m_cur->vert));
			Eigen::Vector3d n2(editMesh.get_vnormal(editMesh.next(*iter.m_cur).vert));
			Eigen::Vector3d n3(editMesh.get_vnormal(editMesh.prev(*iter.m_cur).vert));
			v1 += 0.001 * n1;
			v2 += 0.001 * n2;
			v3 += 0.001 * n3;

			point_vecs.push_back(glm::vec3(v1[0], v1[1], v1[2]));
			point_vecs.push_back(glm::vec3(v2[0], v2[1], v2[2]));
			point_vecs.push_back(glm::vec3(v3[0], v3[1], v3[2]));

			isColored[iter.m_cur->face] = true;

		} while(editMesh.advance_iterator(iter));
	}
	
	// data format transformation
	int size = point_vecs.size();
	GLfloat *vertices = new GLfloat[3*size];
	//GLfloat *colours  = new GLfloat[3*size];
	for (int i = 0; i < size; i++)
	{
		vertices[3*i]   = point_vecs[i].x;
		vertices[3*i+1] = point_vecs[i].y;
		vertices[3*i+2] = point_vecs[i].z;
		//colours[3*i]    = colour_vecs[i].x;
		//colours[3*i+1]  = colour_vecs[i].y;
		//colours[3*i+2]  = colour_vecs[i].z;
	}

	// activate and specify pointer to vertex array
	// this could be much more efficient if I wasn't lazy
	GLuint tmp;
	glGenBuffers(1, &tmp);
	glBindBuffer(GL_ARRAY_BUFFER, tmp);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*size, vertices, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	//// draw the lines
	//glLineWidth(3.0);
	//glEnable(GL_LINE_SMOOTH);
	//glDrawArrays(GL_LINES, 0, 3*size);
	//glDisable(GL_LINE_SMOOTH);

	// Draw the triangles
	glDrawArrays(GL_TRIANGLES, 0,  3*size);

	// deactivate vertex arrays after drawing
	glDisableVertexAttribArray(0);
	//glDisableVertexAttribArray(1);
	glDeleteBuffers(1, &tmp);
	//glDeleteBuffers(1, &tmp2);
	delete[] vertices;
	//delete[] colours;
}


void drawVector(int x1, int y1, int x2, int y2)
{
	// vector has zero magnitude
	if (x1 == x2 && y1 == y2)
		return;

	float fx1 = ((float)x1/c_state.width)*2 - 1;
	float fx2 = ((float)x2/c_state.width)*2 - 1;
	float fy1 = ((float)y1/c_state.height)*2 - 1;
	float fy2 = ((float)y2/c_state.height)*2 - 1;

	GLfloat vertices[] = {  fx1, fy1, 0.5f,
							fx2, fy2, 0.5f };

	// activate and specify pointer to vertex array
	// this could be much more efficient if I wasn't lazy
	GLuint tmp;
	glGenBuffers(1, &tmp);
	glBindBuffer(GL_ARRAY_BUFFER, tmp);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6, vertices, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	// draw the lines
	glLineWidth(2.0);
	glEnable(GL_LINE_SMOOTH);
	//w_state->loadColorMaterial(glm::vec4(0, 0, 1, 1));
	glDrawArrays(GL_LINES, 0, 6);
	glDisable(GL_LINE_SMOOTH);

	//// draw the end point as sphere
	//size_t space(10), radius(100);
	//size_t vertexcount((180 / space) * ( 360 / space) * 2);
	//double* vertex = new double [3 * vertexcount];

	//size_t n=0;
	//for(size_t i=0; i <= 180-space; i += space)
	//{
	//	for(size_t j=0; j<=360-space; j += space)
	//	{

	//		vertex[n + 0]=radius*(sin((i*PI)/180))*(sin((j*PI)/180))-fx2;
	//		vertex[n + 1]=radius*(cos((j*PI)/180))*(sin((i*PI)/180))-fy2;
	//		vertex[n + 2]=radius*(cos((i*PI)/180))-0.5f;

	//		vertex[n + 3]=radius*(sin((i*PI)/180))*(sin(((j+space)*PI)/180))-fx2;
	//		vertex[n + 4]=radius*(cos(((j+space)*PI)/180))*(sin((i*PI)/180))-fy2;
	//		vertex[n + 5]=radius*(cos((i*PI)/180))-0.5f;

	//		n += 6;
	//	}
	//}

	//GLuint tmp2;
	//glGenBuffers(1, &tmp2);
	//glBindBuffer(GL_ARRAY_BUFFER, tmp2);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*vertexcount, vertex, GL_DYNAMIC_DRAW);
	//glEnableVertexAttribArray(1);
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//
	////glDrawArrays(GL_TRIANGLE_STRIP, 0, 3 * vertexcount);
	//glBegin(GL_TRIANGLE_STRIP);
	//for(size_t i = 0;i <= vertexcount; i += 3)
	//{
	//	glVertex3f(vertex[i + 0],vertex[i + 1],vertex[i + 2]);
	//}
	//glEnd();


	// deactivate vertex arrays after drawing
	glDisableVertexAttribArray(0);
	glDeleteBuffers(1, &tmp);
	//glDisableVertexAttribArray(1);
	//glDeleteBuffers(1, &tmp2);
	//delete[] vertex;
}


void drawBall(EditMesh& editMesh, size_t vert)
{
	// extract the vertices' position and face normal data
	std::vector<glm::vec3> point_vecs;
	vface_iterator iter; // initial the iterator visiting the one ring of "vertex"
	if (!editMesh.init_iterator(iter, vert))
		std::cerr << "error because of float point" << std::endl;
	do {
		Eigen::Vector3d v1(editMesh.get_vertex(iter.m_cur->vert));
		Eigen::Vector3d v2(editMesh.get_vertex(editMesh.next(*iter.m_cur).vert));
		Eigen::Vector3d v3(editMesh.get_vertex(editMesh.prev(*iter.m_cur).vert));

		// offset the vertices a bit
		Eigen::Vector3d n1(editMesh.get_vnormal(iter.m_cur->vert));
		Eigen::Vector3d n2(editMesh.get_vnormal(editMesh.next(*iter.m_cur).vert));
		Eigen::Vector3d n3(editMesh.get_vnormal(editMesh.prev(*iter.m_cur).vert));
		v1 += 0.002 * n1;
		v2 += 0.002 * n2;
		v3 += 0.002 * n3;

		point_vecs.push_back(glm::vec3(v1[0], v1[1], v1[2]));
		point_vecs.push_back(glm::vec3(v2[0], v2[1], v2[2]));
		point_vecs.push_back(glm::vec3(v3[0], v3[1], v3[2]));

	} while(editMesh.advance_iterator(iter));


	// data format transformation
	int size = point_vecs.size();
	GLfloat *vertices = new GLfloat[3*size];
	//GLfloat *colours  = new GLfloat[3*size];
	for (int i = 0; i < size; i++)
	{
		vertices[3*i]   = point_vecs[i].x;
		vertices[3*i+1] = point_vecs[i].y;
		vertices[3*i+2] = point_vecs[i].z;
	}

	// activate and specify pointer to vertex array
	// this could be much more efficient if I wasn't lazy
	GLuint tmp;
	glGenBuffers(1, &tmp);
	glBindBuffer(GL_ARRAY_BUFFER, tmp);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*size, vertices, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	//// draw the lines
	//glLineWidth(3.0);
	//glEnable(GL_LINE_SMOOTH);
	//glDrawArrays(GL_LINES, 0, 3*size);
	//glDisable(GL_LINE_SMOOTH);

	// Draw the triangles
	glDrawArrays(GL_TRIANGLES, 0,  3*size);

	glDisableVertexAttribArray(0);
	glDeleteBuffers(1, &tmp);
	delete[] vertices;
	
}

#endif // MESH_UTILS_H