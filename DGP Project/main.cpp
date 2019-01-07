/* Copyright (c) Russell Gillette
 * December 2013
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell co
 es of the Software, and
 * to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define GLFW_INCLUDE_GLU

#ifdef _WIN32
#  include "GL/glew.h"
#  include "GLFW/glfw3.h"
# elif __APPLE__
#  include <GL/glew.h>
#  include <GLFW/glfw3.h>
#else
#  include <GL/glew.h>
#  include <GLFW/glfw3.h>
#endif

#include <iostream>

#include <stdio.h>
#include <math.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>

#include "ShaderUtils.h"
#include "ControlState.h"
#include "WorldState.h"
#include "RenderState.h"
#include "DrawMesh.h"
#include "MeshUtils.h"
#include "TextureUtils.h"
#include "TextureTypes.h"
#include "Utils.h"

#include "EditMesh.h"

WorldState *w_state;
RenderState *r_state[2];
// NOTE: defined in ControlState.h 
// ControlState c_state;

int mesh_file_size = 10;
int mesh_curr = 0;
char *mesh_files[] = {//"Mesh/cat.obj", "Mesh/pyramid.obj", "Mesh/rabbit.obj",
						"Mesh/cube.obj", "Mesh/sphere.obj", "Mesh/camel_simple.obj", "Mesh/camel.obj","Mesh/cow_head.obj",
						"Mesh/cow1.obj", "Mesh/horse.obj" , "Mesh/cow2.obj", "Mesh/octopus.obj", "Mesh/icosahedron.obj"};


Model    *g_mesh;
DrawMesh *g_axis; // NOTE: only a single axis
DrawMesh *g_prioVert;
bool isColorPriorVert(true);
bool isColorDeformVert(true);
bool isRegionSelect(false);
bool isControlSelect(false);
bool isControlMove(false);
bool isSelectThrough(false);



void get_select_vert();
void get_select_region();
void get_move_vert();
void state_clear();
void drawMeshTemp();
void drawControlVertex();
void drawInfluenceRegion();


// the display loop, where all of the code that actually
// changes what you see goes
void display()
{
    /* limit framerate to 60 fps */
    double curr = 0;
    if ((curr = glfwGetTime()) < 0.016666667) // curr < ~ 1/60
        return;

    // start counting over
    glfwSetTime(0.0);

    // Clear the buffer we will draw into.
    glClearColor(0.549, 0.47, 0.937, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Setup camera
    w_state->resetProjection(c_state.aspectRatio());
    w_state->updateView(c_state.viewTheta, c_state.viewPhi,
                        c_state.viewDepth, c_state.viewPan);

    /***********************************
     * Apply pending filter
     ***********************************/
    /* CS524_INPUT 
     * you can add key bindings for your functions through here and c_state */
    switch(c_state.op)
    {
		case EDIT_BUTTERY_SUBDIV :	// for butterfly subdivision on closed mesh
			{
				if (!g_mesh->edit_subdivision())
					std::cerr << "error in subdividing!" << std::endl;	
				c_state.op = EDIT_NONE;
			} break;
		case EDIT_Loop_SUBDIV :	// for loop subdivision on closed mesh
			{
				if (!g_mesh->edit_subdivision(1))
					std::cerr << "error in subdividing!" << std::endl;	
				c_state.op = EDIT_NONE;
			} break;
		
		case EDIT_SIMPLIFY_NUM :	// for mesh simplification, setting number
			{
				// error select
				std::cout << "Please select error metric: 0 for point-to-plane distance, 1 for 1-ring angle" << endl;
				std::size_t mode(1);
				std::cin >> mode;
				if (mode != 0 && mode != 1)
					std::cout << "0 for point-to-plane distance, 1 for 1-ring angle" << endl;

				if (!g_mesh->edit_set_error_simplify(mode))
					std::cerr << "error in simplifying!" << std::endl;	
				
				// number set
				std::cout << "Please key in the number of vertices to be removed" << endl;
				std::size_t num(1);
				std::cin >> num;
				std::cout << "Key right for simplify forward once and left for backward" << endl;
				if (!g_mesh->edit_set_num_simplify(num))
						std::cerr << "error in simplifying!" << std::endl;	
				c_state.op = EDIT_NONE;
			} break;

		case EDIT_SIMPLIFY_ERROR :	// for mesh simplification, setting number
			{

				std::cout << "Please select error metric: 1 for 1-ring angle, 0 for point-to-plane distance" << endl;
				std::size_t num(1);
				std::cin >> num;
				if (num != 0 && num != 1)
					std::cout << "1 for 1-ring angle, 2 for point-to-plane distance" << endl;

				if (!g_mesh->edit_set_error_simplify(num))
					std::cerr << "error in simplifying!" << std::endl;	
				c_state.op = EDIT_NONE;
			} break;

		case EDIT_SIMPLIFY_COLORING :	// for mesh simplification, setting number
			{
				if (isColorPriorVert)
					isColorPriorVert = false;
				else
					isColorPriorVert = true;
				
				c_state.op = EDIT_NONE;
			} break;


		//--------DEFORMATION CODE------------//
		case EDIT_DEFORM_COLORING :	// coloring
			{
				if (isColorDeformVert)
					isColorDeformVert = false;
				else
					isColorDeformVert = true;

				c_state.op = EDIT_NONE;
			} break;


		case EDIT_DEFORM_VERT_SELECT :	// select the control vert 
			{
				if (!isControlSelect)
				{
					isControlSelect = true;
					g_mesh->get_editMesh()->clear_control_vert();
					isRegionSelect = false;
					isControlMove = false;
					std::cout << "Select the control vertex." << endl;
					
				}
				else
				{
					isControlSelect = false;
					std::cout << "control vertex selection completed." << endl;
				}	

				c_state.op = EDIT_NONE;
			} break;

		case EDIT_DEFORM_ROI_SELECT :	// select the ROI
			{
				if (!isRegionSelect)
				{
					isRegionSelect = true;
					g_mesh->get_editMesh()->clear_regionOfInf_vert();
					if (isControlSelect)
					isControlSelect = false;
					isControlMove = false;
					std::cout << "select ROI" << endl;

					
				}
				else
					{
						isRegionSelect = false;
						std::vector<size_t> selected_reg_influence_verts;
						g_mesh->get_editMesh()->get_regionOfInf_vert(selected_reg_influence_verts);
						if (selected_reg_influence_verts.empty())
						std::cout << "No ROI was selected!" << endl;
						else
						std::cout << "ROI selection completed." << endl;
					}

				c_state.op = EDIT_NONE;
			} break;

		case EDIT_DEFORM_VERT_MOVE :	// drag the control vert
			{
				if (isControlMove)
				{
					isControlMove = false;
					size_t control_vert;
					g_mesh->get_editMesh()->get_control_vert(control_vert);
					if (control_vert == std::numeric_limits<size_t>::max())
						std::cout << "No control vertex was selected!" << endl;
					else
						std::cout << "drag the control vertex to start the shape deformation" << endl;
				}
				else
				{
					if (isRegionSelect || isControlSelect)
						std::cout << "You need to select the ROI firt and then you can select the control vertex!" << endl;
					else
						{
							isControlMove = true;
							g_mesh->get_editMesh()->clear_displacementQuanta();
							std::cout << "Drag the control vertex!" << endl;
						}
				}	

				c_state.op = EDIT_NONE;
			} break;

		case EDIT_DEFORM_ROTATION:	// include rotation in the algorithm
			{
				//state_clear();
				bool& IncludeRot = g_mesh->get_editMesh()->IncludeRot;
				if (IncludeRot)
					IncludeRot = false;
				else 
					IncludeRot = true;

				c_state.op = EDIT_NONE;
			} break;

		case EDIT_DEFORM_REGION_SELECT_THROUGH:	// 
			{
				if (isSelectThrough)
				{
					isSelectThrough = false;
				}
				else
				{
					isSelectThrough = true;
				}
				c_state.op = EDIT_NONE;
			} break;

		case EDIT_REMESH:
			{

				std::cout <<"start remeshing"<<endl;
				g_mesh->edit_remesh();
				c_state.op = EDIT_NONE;
			} break;
		case EDIT_LOADMESH:
			{

				std::cout <<"perfomr marching cubes "<<endl;
				g_mesh->edit_loadMesh();
				c_state.op = EDIT_NONE;
			} break;
		case EDIT_LOAD30:
			{

				std::cout <<"load the mesh generated by marshing cubes- 30 slices"<<endl;
				g_mesh->edit_load30Slices();
				c_state.op = EDIT_NONE;
			} break;
		case EDIT_LOAD56:
			{

				std::cout <<"load the mesh generated by marshing cubes- 56 slices"<<endl;
				g_mesh->edit_load56Slices();
				c_state.op = EDIT_NONE;
			} break;

		case EDIT_LOAD113:
			{

				std::cout <<"load the mesh generated by marshing cubes- 113 slices"<<endl;
				g_mesh->edit_load113Slices();
				c_state.op = EDIT_NONE;
			} break;

		case EDIT_NONE:
		default: 
			break;
    }

	// control the simplification
	if (c_state.arrR)
	{
		g_mesh->edit_forward_simplify();
		c_state.arrR = false;
	}

	if (c_state.arrL)
	{
		g_mesh->edit_backward_simplify();
		c_state.arrL = false;
	}

	// Load meshes in a loop
    if (c_state.reload)
    {
        state_clear();
		delete g_mesh;
        g_mesh = loadModelFromFile(*r_state[0], mesh_files[mesh_curr]);
        mesh_curr = (mesh_curr + 1) % mesh_file_size;
        c_state.reload = false;
    }

    /***********************************
     * XYZ Axis Code
     ***********************************/
    w_state->useProgram(0);
    w_state->loadLights();

    //Draw X axis in red
    w_state->loadColorMaterial(glm::vec4(1, 0, 0, 1));
    w_state->loadObjectTransforms(glm::rotate(glm::mat4(),-90.0f, glm::vec3(0, 0, 1)));
    g_axis->drawMesh();

    //Draw Y axis in green
    w_state->loadColorMaterial(glm::vec4(0, 1, 0, 1));
    w_state->loadTransforms();
    g_axis->drawMesh();

    //Draw Z axis in blue
    w_state->loadColorMaterial(glm::vec4(0, 0, 1, 1));
    w_state->loadObjectTransforms(glm::rotate(glm::mat4(),90.0f, glm::vec3(1, 0, 0)));
    g_axis->drawMesh();


	if (!isControlMove)
	{
		/***********************************
		    * Mesh Code
		 ***********************************/
		w_state->useProgram(1);
		// load values into shader
		glm::vec2 screen(c_state.width, c_state.height);
		
		glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
		glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
		w_state->loadTransforms();
		w_state->loadMaterials();
		w_state->loadLights();
		w_state->loadTextures();
		g_mesh->drawMesh();
	}


	////////////////////////////////////////////////
	//// draw priority vertices
	w_state->useProgram(0);
	w_state->loadLights();

	// initialize the  g_prioVert
	std::vector<size_t> vertVecColor;
	g_mesh->get_editMesh()->get_color_vert(vertVecColor);
	EditMesh& editMesh(*g_mesh->get_editMesh());
	if (!vertVecColor.empty() && isColorPriorVert)
	{
		std::vector<size_t> vertVecFront;
		vertVecFront.push_back(vertVecColor.front());
		vertVecColor.erase(vertVecColor.begin());
	
		glm::vec2 screen(c_state.width, c_state.height);
		glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
		glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
		w_state->loadTransforms();
		w_state->loadLights();
		w_state->loadColorMaterial(glm::vec4(1, 0, 0, 1));
		drawTriangles(*g_mesh->get_editMesh(), vertVecFront);

		w_state->loadColorMaterial(glm::vec4(0, 0, 1, 1));
		drawTriangles(*g_mesh->get_editMesh(), vertVecColor);
		//drawLines(_points, _colour_vecs);
	}

	/////////////////////////////////////////////////////////////////
	//draw deformation


	// get influence region
	if (isRegionSelect)
		get_select_region();

	if (isControlSelect)
		get_select_vert();

	if (isControlMove)
		get_move_vert();

	if (isColorDeformVert)
	{
		drawControlVertex();
		drawInfluenceRegion();
	}


	//// draw selected verts
	//w_state->useProgram(0);
	//w_state->loadLights();

	//// draw influence region
	//std::vector<size_t> selected_reg_influence_verts;
	//g_mesh->get_editMesh()->get_region_vert(selected_reg_influence_verts);
	////EditMesh& editMesh(*g_mesh->get_editMesh());
	//if (!selected_reg_influence_verts.empty() && isColorDeformVert)
	//{
	//	glm::vec2 screen(c_state.width, c_state.height);
	//	glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
	//	glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
	//	w_state->loadTransforms();
	//	w_state->loadLights();
	//	w_state->loadColorMaterial(glm::vec4(0.0, 0.7, 0.3, 0.1));
	//	drawTriangles(*g_mesh->get_editMesh(), selected_reg_influence_verts);

	//	drawMeshTemp();
	//}

	//// draw control vertex
	//size_t control_vert;
	//g_mesh->get_editMesh()->get_control_vert(control_vert);
	//if (control_vert != std::numeric_limits<size_t>::max() && isColorDeformVert)
	//{
	//	glm::vec2 screen(c_state.width, c_state.height);
	//	glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
	//	glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
	//	w_state->loadTransforms();
	//	w_state->loadLights();
	//	w_state->loadColorMaterial(glm::vec4(1.0, 0.0, 0.0, 0.1));
	//	drawBall(*g_mesh->get_editMesh(), control_vert);

	//	if(!g_mesh->get_editMesh()->boundaryColor.empty())
	//	{
	//		for (size_t i(0); i < g_mesh->get_editMesh()->boundaryColor.size(); ++i)
	//		{
	//			drawBall(*g_mesh->get_editMesh(), g_mesh->get_editMesh()->boundaryColor[i]);
	//		}
	//	}

	//	//drawMeshTemp();
	//}


	//////////////////////////////////////////////
 //   
 //   /*************************************
 //    * Draw Selection Box
 //    *************************************/
 //   w_state->useProgram(2);
 //   glm::vec3 s_bl, s_tr;
 //   c_state.getMouseSelection(s_bl, s_tr);

 //   glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"bot_left"), 1, glm::value_ptr(s_bl));
 //   glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"top_right"), 1, glm::value_ptr(s_tr));

 ////   drawBox(s_bl.x, s_bl.y, s_tr.x, s_tr.y);

    glfwSwapBuffers(c_state.window);
    glfwPollEvents();
}

void drawMeshTemp()
{
	
	/***********************************
	    * Mesh Code
	 ***********************************/
	w_state->useProgram(1);
	// load values into shader
	glm::vec2 screen(c_state.width, c_state.height);
	
	glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
	glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
	w_state->loadTransforms();
	w_state->loadMaterials();
	w_state->loadLights();
	w_state->loadTextures();
	g_mesh->drawMesh();
}

void drawInfluenceRegion()
{
	// draw selected verts
	w_state->useProgram(0);
	w_state->loadLights();

	// draw influence region
	std::vector<size_t> selected_reg_influence_verts;
	g_mesh->get_editMesh()->get_regionOfInf_vert(selected_reg_influence_verts);
	//EditMesh& editMesh(*g_mesh->get_editMesh());
	if (!selected_reg_influence_verts.empty() && isColorDeformVert)
	{
		glm::vec2 screen(c_state.width, c_state.height);
		glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
		glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
		w_state->loadTransforms();
		w_state->loadLights();
		w_state->loadColorMaterial(glm::vec4(1.0, 0.0, 0.0, 0.1));
		drawTriangles(*g_mesh->get_editMesh(), selected_reg_influence_verts);

		//drawMeshTemp();
	}
}

void drawControlVertex()
{

	// draw selected verts
	w_state->useProgram(0);
	w_state->loadLights();

	// draw control vertex
	size_t control_vert;
	g_mesh->get_editMesh()->get_control_vert(control_vert);
	if (control_vert != std::numeric_limits<size_t>::max() && isColorDeformVert)
	{
		glm::vec2 screen(c_state.width, c_state.height);
		glUniform1i(glGetUniformLocation(w_state->getCurrentProgram(), "view_mode"), c_state.view_mode);
		glUniform2fv(glGetUniformLocation(w_state->getCurrentProgram(), "scale"), 1, glm::value_ptr(screen));
		w_state->loadTransforms();
		w_state->loadLights();
		w_state->loadColorMaterial(glm::vec4(0.0, 0.5, 0.5, 0.1));
		drawBall(*g_mesh->get_editMesh(), control_vert);

		//if(!g_mesh->get_editMesh()->boundaryColor.empty())
		//{
		//	for (size_t i(0); i < g_mesh->get_editMesh()->boundaryColor.size(); ++i)
		//	{
		//		drawBall(*g_mesh->get_editMesh(), g_mesh->get_editMesh()->boundaryColor[i]);
		//	}
		//}

		//drawMeshTemp();
	}
}

// move the control vertex to the closest point on the project ray
void get_move_vert()
{
	w_state->useProgram(2);
	glm::vec3 s_bl(0, 0, 0.5f), s_tr(0, 0, 0.5f);
	c_state.getMouseStartEnd(s_bl, s_tr);

	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"bot_left"), 1, glm::value_ptr(s_bl));
	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"top_right"), 1, glm::value_ptr(s_tr));

	if (s_bl != glm::vec3(0, 0, 0.5f) || s_tr != glm::vec3(0, 0, 0.5f))
	{
		drawVector(s_bl.x, s_bl.y, s_tr.x, s_tr.y);


		// get the displacement vector from editMesh
		Eigen::Vector3d& DisplacementQuanta(g_mesh->get_editMesh()->DisplacementQuanta);		

		// get the normal of the control vertex
		size_t& control_vert(g_mesh->get_editMesh()->control_vert);
		if (control_vert == std::numeric_limits<size_t>::max())
		{
				std::cout << "Please select control vertex first!" << endl;
				return;
		}

		// project the control and normal to screen
		Eigen::Vector3d vert(g_mesh->info_vertex(control_vert));

		// get projection ray and project point of the end position of control vert
		GLint select_viewport[4];
		glGetIntegerv(GL_VIEWPORT,select_viewport);
		glm::vec3 endPoint = glm::unProject(glm::vec3(s_tr.x,s_tr.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 endRay   = glm::unProject(glm::vec3(s_tr.x,s_tr.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));

		// compute the closest point on the ray
		Eigen::Vector3d ray(endRay.x, endRay.y, endRay.z);
		Eigen::Vector3d p0(endPoint.x, endPoint.y, endPoint.z);
		ray.normalize();

		DisplacementQuanta = -(vert - p0 - ray.dot(vert - p0) * ray);

		if (DisplacementQuanta[0] == std::numeric_limits<double>::max())
		{
			std::cout << "Error occurred when computing the displacement of control vertex!" << endl;
			return;
		}

		g_mesh->get_editMesh()->move_control_vert();
		drawMeshTemp();

		if (isColorDeformVert)
		{
			drawControlVertex();
			drawInfluenceRegion();
		}
		
	}
	else
	{
		drawMeshTemp();

		if (isColorDeformVert)
		{
			drawControlVertex();
			drawInfluenceRegion();
		}
	}

}

void get_select_vert()
{
	w_state->useProgram(2);
	glm::vec3 s_bl, s_tr;
	c_state.getMouseSelection(s_bl, s_tr);

	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"bot_left"), 1, glm::value_ptr(s_bl));
	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"top_right"), 1, glm::value_ptr(s_tr));

	drawBox(s_bl.x, s_bl.y, s_tr.x, s_tr.y);

	if (s_bl != glm::vec3(0, 0, 0) || s_tr != glm::vec3(0, 0, 0))
	{
		// get the container to store region verts
		size_t& control_vert(g_mesh->get_editMesh()->control_vert);
		//Eigen::Vector3d normal(g_mesh->info_vertexNormal(control_vert));

		//determine which vertices are in the selection box
		GLint select_viewport[4];
		glGetIntegerv(GL_VIEWPORT,select_viewport);
		glm::vec3 bl     = glm::unProject(glm::vec3(s_bl.x,s_bl.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 bl_ray = glm::unProject(glm::vec3(s_bl.x,s_bl.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 br     = glm::unProject(glm::vec3(s_tr.x,s_bl.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 br_ray = glm::unProject(glm::vec3(s_tr.x,s_bl.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tr     = glm::unProject(glm::vec3(s_tr.x,s_tr.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tr_ray = glm::unProject(glm::vec3(s_tr.x,s_tr.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tl     = glm::unProject(glm::vec3(s_bl.x,s_tr.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tl_ray = glm::unProject(glm::vec3(s_bl.x,s_tr.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));

		int vert_size = g_mesh->info_sizev();
		for (int i = 0; i < vert_size; i++)
		{
			Eigen::Vector3d vert = g_mesh->info_vertex(i);
			Eigen::Vector3d normal(g_mesh->info_vertexNormal(i));

			if (vert_inside_select_box(Eigen::Vector3d(bl.x,bl.y,bl.z),
				Eigen::Vector3d(bl_ray.x,bl_ray.y,bl_ray.z),
				Eigen::Vector3d(br.x,br.y,br.z),
				Eigen::Vector3d(br_ray.x,br_ray.y,br_ray.z),
				Eigen::Vector3d(tr.x,tr.y,tr.z),
				Eigen::Vector3d(tr_ray.x,tr_ray.y,tr_ray.z),
				Eigen::Vector3d(tl.x,tl.y,tl.z),
				Eigen::Vector3d(tl_ray.x,tl_ray.y,tl_ray.z),
				vert))
			{
				// decide the side: front or back of vertex control_vert
				Eigen::Vector3d n1(bl_ray.x, bl_ray.y, bl_ray.z);
				Eigen::Vector3d n2(br_ray.x, br_ray.y, br_ray.z);
				Eigen::Vector3d n3(tr_ray.x, tr_ray.y, tr_ray.z);
				Eigen::Vector3d n4(tl_ray.x, tl_ray.y, tl_ray.z);

				n1.normalize();
				n2.normalize();
				n3.normalize();
				n4.normalize();

				if (n1.dot(n2) < 0 || n1.dot(n3) < 0 || n1.dot(n4) < 0 || n2.dot(n3) < 0 || n2.dot(n3) < 0 || n3.dot(n4) < 0)
				{
					std::cout << "The computed ray do not face the same side!" << endl;
				}

				Eigen::Vector3d n0 = (n1 + n2 + n3 + n4) / 4.0;
				n0.normalize();

				if (n0.dot(normal) > 0)
					continue;
				else
					control_vert = i;
			}
		}

	}

	drawControlVertex();
}

void get_select_region()
{
	w_state->useProgram(2);
	glm::vec3 s_bl, s_tr;
	c_state.getMouseSelection(s_bl, s_tr);
	//c_state.getSelectSquare(s_bl, s_tr);

	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"bot_left"), 1, glm::value_ptr(s_bl));
	glUniform3fv(glGetUniformLocation(w_state->getCurrentProgram(),"top_right"), 1, glm::value_ptr(s_tr));

	drawBox(s_bl.x, s_bl.y, s_tr.x, s_tr.y);

	if (s_bl != glm::vec3(0, 0, 0) || s_tr != glm::vec3(0, 0, 0))
	{
		// get the container to store region verts
		std::set<size_t>& selected_reg_influence_verts(g_mesh->get_editMesh()->selected_reg_influence_verts);

		//determine which vertices are in the selection box
		GLint select_viewport[4];
		glGetIntegerv(GL_VIEWPORT,select_viewport);
		glm::vec3 bl     = glm::unProject(glm::vec3(s_bl.x,s_bl.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 bl_ray = glm::unProject(glm::vec3(s_bl.x,s_bl.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 br     = glm::unProject(glm::vec3(s_tr.x,s_bl.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 br_ray = glm::unProject(glm::vec3(s_tr.x,s_bl.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tr     = glm::unProject(glm::vec3(s_tr.x,s_tr.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tr_ray = glm::unProject(glm::vec3(s_tr.x,s_tr.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tl     = glm::unProject(glm::vec3(s_bl.x,s_tr.y,0), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));
		glm::vec3 tl_ray = glm::unProject(glm::vec3(s_bl.x,s_tr.y,1), w_state->view * w_state->model, w_state->projection, glm::vec4(0.0,0.0,select_viewport[2],select_viewport[3]));

		int vert_size = g_mesh->info_sizev();
		for (int i = 0; i < vert_size; i++)
		{
			Eigen::Vector3d vert = g_mesh->info_vertex(i);
			Eigen::Vector3d normal(g_mesh->info_vertexNormal(i));
			if (vert_inside_select_box(Eigen::Vector3d(bl.x,bl.y,bl.z),
				Eigen::Vector3d(bl_ray.x,bl_ray.y,bl_ray.z),
				Eigen::Vector3d(br.x,br.y,br.z),
				Eigen::Vector3d(br_ray.x,br_ray.y,br_ray.z),
				Eigen::Vector3d(tr.x,tr.y,tr.z),
				Eigen::Vector3d(tr_ray.x,tr_ray.y,tr_ray.z),
				Eigen::Vector3d(tl.x,tl.y,tl.z),
				Eigen::Vector3d(tl_ray.x,tl_ray.y,tl_ray.z),
				vert))
			{
				if (isSelectThrough)
				{
					selected_reg_influence_verts.insert(i);
				}
				else
				{
					// decide the side: front or back of vertex control_vert
					Eigen::Vector3d n1(bl_ray.x, bl_ray.y, bl_ray.z);
					Eigen::Vector3d n2(br_ray.x, br_ray.y, br_ray.z);
					Eigen::Vector3d n3(tr_ray.x, tr_ray.y, tr_ray.z);
					Eigen::Vector3d n4(tl_ray.x, tl_ray.y, tl_ray.z);

					n1.normalize();
					n2.normalize();
					n3.normalize();
					n4.normalize();

					if (n1.dot(n2) < 0 || n1.dot(n3) < 0 || n1.dot(n4) < 0 || n2.dot(n3) < 0 || n2.dot(n3) < 0 || n3.dot(n4) < 0)
					{
						std::cout << "The computed ray do not face the same side!" << endl;
					}

					Eigen::Vector3d n0 = (n1 + n2 + n3 + n4) / 4.0;
					n0.normalize();

					if (n0.dot(normal) > 0)
						continue;
					else
						selected_reg_influence_verts.insert(i);

				}
				
			}
		}

	}

	drawInfluenceRegion();
}

void state_clear()
{
	//isColorPriorVert = true;
	//isColorDeformVert = true;
	isRegionSelect = false;
	isControlSelect = false;
	isControlMove = false;
	isSelectThrough = false;

	g_mesh->get_editMesh()->clear_control_vert();
	g_mesh->get_editMesh()->clear_regionOfInf_vert();
	g_mesh->get_editMesh()->clear_displacementQuanta();
}

// setup
int main(int argc, char *argv[])
{
	// Test
	//EditMesh::test();

    GLenum err = 0;
    /*********************************************
     * GLFW SETUP
     *********************************************/
    err = glfwInit();
    if (!err)
    {
        fputs("Failed to load the GLFW library", stderr);
        exit(EXIT_FAILURE);
    }

#ifdef __APPLE__
     glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	 glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	 glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	 glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif
    
    /*********************************************
     * STATE SETUP (initialize gl context)
     *********************************************/
    // must be setup before glew so that a valid openGL
    // context exists (created with the window)

    w_state = new WorldState();
    c_state.init(*w_state);

    /*********************************************
     * GLEW SETUP
     *********************************************/
#ifdef __APPLE__
     glewExperimental = GL_TRUE;
#endif
    err = glewInit();
    if (err != GLEW_OK)
    {
        fputs("Failed to initialize the GLEW library", stderr);
        exit(EXIT_FAILURE);
    }

    /*********************************************
     * STATE SETUP (construct render states)
     *********************************************/
    // must be setup after glew so that GL array
    // objects exist

    r_state[0] = new RenderState();
    r_state[1] = new RenderState();

    /*********************************************
     * SHADER SETUP
     *********************************************/
    // read default shaders from file
    GLuint shaderProgram[3] = {0};
    GLuint shaders[3] = {0};

    buildShader(GL_VERTEX_SHADER, "axes.vs.glsl", shaders[0]);
    buildShader(GL_FRAGMENT_SHADER, "default.fs.glsl", shaders[1]);

    // create axis shader program
    shaderProgram[0] = buildProgram(2, shaders);

    // create the shaders for the mesh
    buildShader(GL_VERTEX_SHADER,   "mesh.vs.glsl", shaders[0]);
    buildShader(GL_GEOMETRY_SHADER, "wireframe.gs.glsl", shaders[2]);
    buildShader(GL_FRAGMENT_SHADER, "wireframe.fs.glsl", shaders[1]);
    shaderProgram[1] = buildProgram(3, shaders);

    // load shaders to render selection
    buildShader(GL_VERTEX_SHADER,   "passthrough.vs.glsl", shaders[0]);
    buildShader(GL_FRAGMENT_SHADER, "select.fs.glsl", shaders[1]);
    shaderProgram[2] = buildProgram(2, shaders);

    // bind shader program
    w_state->setProgram(0, shaderProgram[0]);
    w_state->setProgram(1, shaderProgram[1]);
    w_state->setProgram(2, shaderProgram[2]);
    w_state->useProgram(0);

    // setup the transform matrices and uniform variables
    w_state->loadTransforms();
    w_state->loadLights();
    w_state->loadMaterials();


	/*********************************************
     * LOAD MESH
     *********************************************/
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/cube.obj");
	//g_mesh = loadModelFromFile(*r_state[0],"Mesh/pyramid.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/icosahedron.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/sphere.obj");
	
	//g_mesh = loadModelFromFile(*r_state[0], "debug.obj");
    //g_mesh = loadModelFromFile(*r_state[0], "Mesh/camel.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/cow1.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/cow_head.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/cow_ear.obj");
    g_mesh = loadModelFromFile(*r_state[0], "Mesh/camel_simple.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/cow2.obj");
	//g_mesh = loadModelFromFile(*r_state[0], "Mesh/octopus.obj");
    //g_mesh = loadModelFromFile(*r_state[0], "Mesh/camel.obj");

    g_axis = createAxis(*r_state[1], 1);

    /*********************************************
     * SET GL STATE
     *********************************************/
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    /*********************************************
     * RENDER LOOP
     *********************************************/
    printHelp();
    glfwSetTime(0.0);
    while (!glfwWindowShouldClose(c_state.window))
        display();

    /*********************************************
     * CLEAN UP
     *********************************************/
    delete g_mesh;
    delete g_axis;

    glfwTerminate();

    exit(EXIT_SUCCESS);
}