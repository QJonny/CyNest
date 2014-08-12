#include "network.h"


// neuron

Neuron::Neuron(int id_, double x_, double y_, double z_):
	id(id_),
	position(x_, y_, z_) {
		connections = vector<int>();
		alpha = ALPHA_THRESHOLD;
		pthread_mutex_init (&mutex , NULL );
		spike_time = 0.0;
		
		// test feature: it randomly selects a neuron
		double n = generateRandomNumber(1, 100);
		if(n >= 10.0 && n < 11.0) {
			selected = true;
		} else {
			selected = false;
		}
}

Neuron::Neuron() {
	Neuron(0, 0.0, 0.0, 0.0);
}

void Neuron::addConnection(int id_) {
	connections.push_back(id_);
}

vector<int>* Neuron::getConnections(){
	return &connections;
}

int Neuron::getId(){
	return id;
}

Vector3d Neuron::getPosition() {
	return position;
}


double Neuron::getAlpha() {
	return alpha;
}

void Neuron::draw() {
	if(id != 0) {
		glColor4f(1.0,1.0,1.0, alpha );
		
		if(alpha > ALPHA_THRESHOLD) { // active
			glColor4f(0.0,1.0,1.0, alpha );
		}
			
		if(selected) {
			glColor4f(1.0,0.0,0.0, alpha );
		} 

		glVertex3f(position.x(), position.y(), position.z());
	}
}



void Neuron::fire(double time_) {
	pthread_mutex_lock (&mutex);
	spikes_buffer.push_back(time_);
	pthread_mutex_unlock (&mutex);
}

void Neuron::update(double time_) {
	if(id != 0) {
		double item = -1.0;
		pthread_mutex_lock (&mutex);
		
		alpha = -(time_ - spike_time) + 1.0;
		
		if(alpha < ALPHA_THRESHOLD) {
			alpha = ALPHA_THRESHOLD;
		}

		// we loop over the spikes list and delete those which are 
		// already expired (if any)
		while(spikes_buffer.size() > 0 && item < time_) {	
			item = spikes_buffer.at(0);
			
			spikes_buffer.erase (spikes_buffer.begin());
		}
		
		if(item != -1.0) {
			// if the spike time is during this simulation interval we fire,
			if(item >= time_ && item < time_ + SIMULATION_DELTA) {
				spike_time = item;
				alpha = 1.0; // fire
			}	
			else { // otherwise we put it inside the list again
				spikes_buffer.insert (spikes_buffer.begin(), item);
			}
		}
		pthread_mutex_unlock (&mutex);
	}
}


void Neuron::select() {
	selected = true;
}

void Neuron::unselect() {
	selected = false;
}

bool Neuron::isSelected() {
	return selected;
}



Vector3d Neuron::get_2D_pos_from_3D(Vector3d pos_, double w, double h) {
	GLdouble modelview[16];
	GLdouble projection[16];
	GLint viewport[16];
	GLdouble x_;
	GLdouble y_;
	GLdouble z_;
	
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    gluProject(pos_.x(),pos_.y(),pos_.z(), modelview, projection, viewport, &x_, &y_, &z_);
    
    return Vector3d(-1.0 + 2.0*x_/w, -1.0 + 2.0*y_/h, z_);
}


bool Neuron::isMouseFocused(double x, double y, double zone_width, double zone_height, double window_width, double window_height) {
	Vector3d pos2d = get_2D_pos_from_3D(position, window_width, window_height);
	
	if(x >= pos2d.x() && x <= pos2d.x() + zone_width && y >= pos2d.y() && y <= pos2d.y() + zone_height) {
		return true;
	}
	
	return false;
}



// camera

void Camera::init() {
	// display parameters initialization
	theta   = 0.0;
	phi     = 0.0;
	pos.set(0.0, 60.0, 0.0);
	lookAtVector.set(0.0, 0.0, 0.0);
	mode = MODE_CENTERED;
}

void Camera::update() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	switch(mode) {
	case MODE_FREE:
		lookAtVector = pos.add(Vector3d(cos(theta) * cos(phi),   sin(theta) * cos(phi),   sin(phi)));
		break;
	case MODE_CENTERED:
		pos.set(pos.norm() * cos(theta) * cos(phi),   pos.norm() * sin(theta) * cos(phi),   pos.norm() * sin(phi));
		break;
	}

	gluLookAt(pos.x(), pos.y(), pos.z(),  lookAtVector.x(), lookAtVector.y(), lookAtVector.z(),  0.0, 0.0, 1.0);
}

void Camera::setMode(int mode_) {
	Vector3d dir = lookAtVector.sub(pos).normalize();
	
	switch(mode_) {
	case MODE_FREE:
		phi = asin(dir.z() / dir.norm());
		theta = atan(dir.y() / dir.x());
		
		// we have to take care of the case in which x<0
		// and y<0, which gives theta > 0 (wrong)
		if(dir.x() < 0) {
			theta = fmod(theta + M_PI, 2*M_PI);
		}
		
		break;
	case MODE_CENTERED:
		phi = asin(pos.z() / pos.norm());
		theta = atan(pos.y() / pos.x());
		
		// we have to take care of the case in which x<0
		// and y<0, which give theta > 0 (wrong)
		if(pos.x() < 0) {
			theta = fmod(theta + M_PI, 2*M_PI);
		}
		
		lookAtVector.set(0.0, 0.0, 0.0);
		break;
	}
	
	mode = mode_;
}

void Camera::up() {
	if(phi < (M_PI / 2) - ANGLE_DIFF) {
		phi = phi + ANGLE_DIFF;
	}
}

void Camera::down() {
	if(phi > -(M_PI / 2) + ANGLE_DIFF) {
		phi = phi - ANGLE_DIFF;
	}
}

void Camera::right() {
	switch(mode) {
	case MODE_FREE:
		theta = fmod(theta - ANGLE_DIFF, 2*M_PI);
		break;
	case MODE_CENTERED:
		theta = fmod(theta + ANGLE_DIFF, 2*M_PI);
		break;
	}
}

void Camera::left() {
	switch(mode) {
	case MODE_FREE:
		theta = fmod(theta + ANGLE_DIFF, 2*M_PI);
		break;
	case MODE_CENTERED:
		theta = fmod(theta - ANGLE_DIFF, 2*M_PI);
		break;
	}
}

void Camera::forward() {
	Vector3d dir = lookAtVector.sub(pos).normalize();
	
	// we add the direction vector to the position
	
	switch(mode) {
	case MODE_FREE:
		pos = pos.add(dir.mul(DIST_DIFF));
		lookAtVector = lookAtVector.add(dir.mul(DIST_DIFF));
		break;
	case MODE_CENTERED:
		if (pos.norm() > DIST_DIFF ) {
			pos = pos.add(dir.mul(DIST_DIFF));
		}
		break;
	}
}

void Camera::backward() {
	Vector3d dir = lookAtVector.sub(pos).normalize();
	
	// we substract the direction vector to the position

	switch(mode) {
	case MODE_FREE:
		pos = pos.sub(dir.mul(DIST_DIFF));
		lookAtVector = lookAtVector.sub(dir.mul(DIST_DIFF));
		break;
	case MODE_CENTERED:
		pos = pos.sub(dir.mul(DIST_DIFF));
		break;
	}
}




