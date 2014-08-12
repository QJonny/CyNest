#include "window.h"

using namespace std;

void Window::init(int width_, int height_, Neuron* neurons_, int* nb_neurons_, int* sim_step, double* sim_time_) {
	neurons = neurons_;
	nb_neurons = nb_neurons_;
	width = width_;
	height = height_;
	simulation_step = sim_step;
	simulation_total_time = sim_time_;


	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER);
	TTF_Init();
	
	surface = SDL_SetVideoMode(width, height, 0, SDL_OPENGL | SDL_DOUBLEBUF );
	    
    init_display();
    
    // background
	glClearColor(0.0, 0.0, 0.0, 1.0);
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_1D);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_SMOOTH);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Window::init_display() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
 
	gluPerspective(45.0, float(width)/float(height), 0.01, 1000.0);
	glMatrixMode(GL_MODELVIEW);
}

void Window::resize(int width_, int height_) {
	width = width_;
	height = height_;
	glViewport(0, 0, width, height);
	
	init_display();
}




void Window::destroy() {
	TTF_Quit();
	SDL_Quit();
}







void Window::update(double time_) {
	double percentage = 100.0;
	if(time_ < *simulation_total_time) {
		percentage = time_ * 100.0 / *simulation_total_time;
	}
	
	sprintf(caption, "CyNEST Graphics Simulator - State: %.1f%, RT Factor: 1/%d",  percentage, *simulation_step);
	SDL_WM_SetCaption(caption, 0);
}


void Window::draw_connections() {
	glPushMatrix();

	glBegin(GL_LINES);
	for(int i=0; i < *nb_neurons; i++) {		
		if(neurons[i].getId() != 0) {
			Vector3d color = Vector3d(1.0, 1.0, 1.0);
			double dAlpha = 0.0;
			
			if(neurons[i].getAlpha() > ALPHA_THRESHOLD) { // actif
				color.set(0.0, 0.5, 0.5);
				dAlpha = 0.4;
			}
			
			if(neurons[i].isSelected()) { // selected
				color.set(1.0, 0.0, 0.0);
				dAlpha = 0.0;
			}

			
			glColor4f(color.x(), color.y(), color.z(), neurons[i].getAlpha() - dAlpha);
			
			Vector3d src = neurons[i].getPosition();
			vector<int>* connections = neurons[i].getConnections();
			for(int j = 0; j < connections->size(); j++) {
				Vector3d dest = neurons[connections->at(j)].getPosition();

				glVertex3f(src.x(),src.y(),src.z());
				glVertex3f(dest.x(),dest.y(),dest.z());
			}
		}
	}
	glEnd();
	glPopMatrix();
}


void Window::draw() {
	// draw points (these should be drawn nicer using shaders)
	glPushMatrix();
	glPointSize(POINT_SIZE);
	glBegin(GL_POINTS);

	for(int i=0; i < *nb_neurons; i++) {
		neurons[i].draw();
	}
	glEnd();
	glPopMatrix();
	
	// draw connections
	draw_connections();
	
	
	glFlush();
	SDL_GL_SwapBuffers();
}
