#ifndef NETWORK_H
#define NETWORK_H


#include "headers.h"

#include "tools.h"


class Neuron
{
private :
	int id;
	Vector3d position;
	double alpha;
	vector<double> spikes_buffer;
	double spike_time;
	pthread_mutex_t mutex;
	bool selected;
	
	vector<int> connections;
	
	Vector3d get_2D_pos_from_3D(Vector3d pos_, double w, double h);
	
public:
	Neuron();
	Neuron(int id_, double x_, double y_, double z_);

	
	void addConnection(int id_);
	
	vector<int>* getConnections();
	int getId();
	double getAlpha();
	Vector3d getPosition();
	
	void fire(double time_);
	
	void update(double time_);
	
	void draw();
	
	void select();
	void unselect();
	bool isSelected();
	
	bool isMouseFocused(double x, double y, double zone_width, double zone_height, double window_width, double window_height);
};


class Camera
{
private:
	double theta;
	double phi;
	
	int mode;
	
	Vector3d lookAtVector;
	Vector3d pos;
public:
	void init();
	void update();
	
	void setMode(int mode_);
	
	void up();
	void down();
	void right();
	void left();
	void forward();
	void backward();
};



#endif


