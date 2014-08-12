#ifndef WINDOW_H
#define WINDOW_H


#include "headers.h"


#include "tools.h"
#include "network.h"





class Window
{
private:
	SDL_Surface *surface;
	Neuron* neurons;
	int* nb_neurons;
	int* simulation_step;
	double* simulation_total_time;
	char caption[200];
	

	
	// Display parameters
	int width;
	int height;

	void init_display();
	
	void draw_connections();

public:
		void init(int width_, int height_, Neuron* neurons_, int* nb_neurons_, int* sim_step, double* sim_time_);
		void destroy();
		
		void resize(int width_, int height_);
		
		void update(double time_);
		
		void draw();
};

#endif
