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

#include "ControlState.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

ControlState c_state = ControlState();

ControlState::~ControlState()
{
    glfwDestroyWindow(window);
}

int ControlState::init(WorldState &w)
{
    this->w = &w;

    width  = 640;
    height = 480;
    /* As of right now we only have one window */
    window = glfwCreateWindow(width, height, "Window", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        fputs("failed to initialize window", stderr);
        return 1; // error
    }
    glfwMakeContextCurrent(window);

    // bind all callbacks
    glfwSetKeyCallback(window, key_callback);
    glfwSetFramebufferSizeCallback(window, reshape_callback);
    glfwSetCursorPosCallback(window, mousePos_callback);
    glfwSetCursorEnterCallback(window, mouseEnter_callback);
    glfwSetMouseButtonCallback(window, mouseBtn_callback);
    glfwSetScrollCallback(window, mouseScroll_callback);

    return 0;
}

int ControlState::deltaArrLR()
{
    return arrR - arrL;
}

int ControlState::deltaArrUD()
{
    return arrD - arrU;
}

void ControlState::clearViewDeltas()
{
    viewTheta = 0;
    viewPhi   = 0;
    viewDepth = 0;
    viewPan   = glm::vec3(0, 0, 0);
}

void printHelp()
{
    printf("==== Help ====\n"
           "___Control___\n"
           "m - switch between modes (view and select)\n\n"
           "___Command___\n"
           "q, esc - exit\n"
           "h - help (you have already figured this out)\n\n"
           "___View___\n"
           "left click and drag      - adjusts view\n"
           "shft+left click and drag - pans the view\n"
           "scroll wheel             - zoom\n");
}

/*****************************************************************************
 * Passive Callback functions
 *****************************************************************************/
// error callback for GLFW
static void error_callback(int error, const char* desc)
{
    fputs(desc, stderr);
}

// callback when window is resized
void reshape_callback(GLFWwindow* window, int w, int h)
{
    c_state.height = h;
    c_state.width  = w;

    glViewport( 0, 0, (GLint)w, (GLint)h );
}

/*****************************************************************************
 * Active Callback functions
 *****************************************************************************/
// callback when a key is pressed
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    /* key code goes here */

    switch(key)
    {
	case GLFW_KEY_LEFT:
		c_state.arrL = (action == GLFW_RELEASE) ? 0 : 1;
		break;
	case GLFW_KEY_RIGHT:
		c_state.arrR = (action == GLFW_RELEASE) ? 0 : 1;
        break;
    case GLFW_KEY_UP:
        c_state.arrU = (action == GLFW_RELEASE) ? 0 : 1;
        break;
    case GLFW_KEY_DOWN:
        c_state.arrD = (action == GLFW_RELEASE) ? 0 : 1;
        break;
    case GLFW_KEY_ESCAPE: // see also Q
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;
	//case GLFW_KEY_D:
	//	if( action == GLFW_RELEASE )
	//		c_state.op = EDIT_DEBUG;
		break;
    case GLFW_KEY_H:
        printHelp();
        break;
    case GLFW_KEY_L:
        c_state.reload = (action == GLFW_RELEASE) ? c_state.reload : 1;
        break;
    case GLFW_KEY_M:
        if (action == GLFW_RELEASE)
            c_state.mode = (RENDER_MODE)((c_state.mode + 1) % MODE_MAX);
        break;
    case GLFW_KEY_N:
        if (action == GLFW_RELEASE) {
            int tmp = (c_state.view_mode + 1);
            c_state.view_mode = (tmp > VIEW_ALL || tmp < VIEW_FACES) ? VIEW_FACES : tmp;
        }
        break;
    case GLFW_KEY_Q:
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;
	///////////////////////////////////////////////////////
    // add by Qiang
	case GLFW_KEY_S:	// for butterfly subdivision on closed mesh
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_BUTTERY_SUBDIV;
		break;
	case GLFW_KEY_U:	// for loop subdivision on closed mesh
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_Loop_SUBDIV;
		break;
	case GLFW_KEY_I:  // for simplification
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_SIMPLIFY_NUM;
		break;

	case GLFW_KEY_O:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_SIMPLIFY_COLORING;
		break;


	// for deformation
	case GLFW_KEY_C:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_COLORING;
		break;

	case GLFW_KEY_X:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_VERT_SELECT;
		break;

	case GLFW_KEY_Z:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_ROI_SELECT;
		break;

	case GLFW_KEY_V:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_VERT_MOVE;
		break;

	case GLFW_KEY_R:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_ROTATION;
		break;

	case GLFW_KEY_T:
		if (action == GLFW_RELEASE)
			c_state.op = EDIT_DEFORM_REGION_SELECT_THROUGH;
		break;
	
	// end for deformation

	case GLFW_KEY_G:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_REMESH;
        break;
	case GLFW_KEY_1:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_LOADMESH;
        break;
	case GLFW_KEY_7:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_LOAD30;
        break;
	case GLFW_KEY_8:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_LOAD56;
        break;
	case GLFW_KEY_9:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_LOAD113;
        break;

	// end of adding

    case GLFW_KEY_0:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_SQRT3_SUBDIV;
        break;
   /* case GLFW_KEY_9:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_RELOCATE_VERTS;
        break;
	case GLFW_KEY_8:
        if (action == GLFW_RELEASE)
            c_state.op = (mode & GLFW_MOD_SHIFT) ? EDIT_COLLAPSE_FAST_EDIT : EDIT_COLLAPSE_EDIT;
        break;
    case GLFW_KEY_7:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_AREA_RELOCATION;
        break;*/
    case GLFW_KEY_6:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_DELAUNAY;
        break;
    case GLFW_KEY_4:
        if (action == GLFW_RELEASE)
            c_state.op = EDIT_PRETTY;
        break;
    }
}

// callback when a mouse button is pressed
static void mouseBtn_callback(GLFWwindow* win, int button, int action, int mod)
{
    /* TODO: any controls relative to pressing the mouse buttons goes here */

    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        c_state.mouseBtnL = (action == GLFW_PRESS) ? 1 : 0;
        c_state.mouseModShft = (mod & GLFW_MOD_SHIFT) ? 1 : 0;
            
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
	{
        c_state.mouseBtnR = (action == GLFW_PRESS) ? 1 : 0;
		if (!c_state.mouseBtnR)
			c_state.select_click = true;
	}
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE)
        c_state.mouseBtnC = (action == GLFW_PRESS) ? 1 : 0;

    if (action == GLFW_PRESS)
        c_state.mouseEvent = true;
    else
    {
       // copy the curser position
		c_state.square_begin = c_state.select_start;
		c_state.square_end = c_state.select_curr;
 
		c_state.select_start = glm::vec3(0, 0, 0);
        c_state.select_curr = glm::vec3(0, 0, 0);		
    }
}

// callback when the mouse is moved. This will be called
// ALOT keep it as light as possible!!
static void mousePos_callback(GLFWwindow* win, double x, double y)
{
    // screen Y coords are inverted.
    y = c_state.height - y;

    // currently used to update camera angles if mouse pressed
    if (c_state.mouseBtnL)
    {
        // Calculate change from last known mouse positon.
        float dx = (float)x - c_state.mouseX;
        float dy = (float)y - c_state.mouseY;
    
        if (c_state.mouseModShft) // update screen pan
        {
            double perc = 50 / c_state.viewDepth;
            c_state.viewPan = c_state.viewPan + glm::vec3(dx/perc, dy/perc, 0);
        }
        else // Update viewing angles.
        {
            c_state.viewTheta = fmod(c_state.viewTheta + 360 + dx / 2, 360);
            c_state.viewPhi   = std::fmin(90.0f, std::fmax(-90.0f, c_state.viewPhi - dy));
        }
    }

    c_state.mouseX = x;
    c_state.mouseY = y;

    if (c_state.mouseBtnR && c_state.mouseEvent)
    {
        c_state.select_start = glm::vec3(x, y, 0);
        c_state.select_curr = c_state.select_start;
		c_state.square_begin = glm::vec3(0, 0, 0);
		c_state.square_end = glm::vec3(0, 0, 0);
        c_state.mouseEvent = false;
    }
    else if (c_state.mouseBtnR)
	{
        c_state.select_curr = glm::vec3(x, y, 0);
		c_state.square_begin = glm::vec3(0, 0, 0);
		c_state.square_end = glm::vec3(0, 0, 0);
		//c_state.select_end = c_state.select_curr;
	}
}

static void mouseScroll_callback(GLFWwindow* win, double x_offset, double y_offset)
{
    // since we would read from mouseScroll, set viewDepth
    // and then clear mouseScroll, I shall refrain from even
    // setting it for this usage
    c_state.viewDepth -= y_offset/30;
}

static void mouseEnter_callback(GLFWwindow* win, int entered)
{
    c_state.mouseInWindow = entered;
}