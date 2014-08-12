#include "graphics_simulator.h"

using namespace std;


// init and end

void GraphicsSimulator::initialize(int port_send, int port_receive, int window_width, int window_height) {
	init_connection(port_send, port_receive);
	window.init(window_width, window_height, neurons, &nb_neurons, &simulation_step, &sim_time);
	camera.init();
	
	simulation_step = INITIAL_SIMULATION_STEP;
	plus_pressed = false;
	minus_pressed = false;
	p_pressed = false;
	k1_pressed = false;
	k2_pressed = false;
}

void GraphicsSimulator::finalize() {
	delete[] neurons;
	listener.destroy();
	sender.destroy();
	window.destroy();
}




// connection management

void GraphicsSimulator::init_connection(int port_send, int port_receive) {
	sender.initiateConnection(port_send, NULL);
	listener.acceptConnection(port_receive, NULL);
	
	sleep(1);
	
	sender.sendMsg("ready", 5);
	char buffer[20];
	double tm[1];
	listener.receiveMsg(buffer, 20);

	if(parseList(buffer, tm)) {	
		sim_time = tm[0];
	}else {
		printf("Problem receiving total simulation time. Program will exit.");
		exit(1);
	}
	sender.sendMsg("ok", 2);
	
	
	receive_positions();
	receive_connections();
}


void GraphicsSimulator::receive_positions() {
	vector<Neuron> neurons_;
	char buffer[50];
	double pos[4];
	int max = 0;

	while(true){
		listener.receiveMsg(buffer, 50); // [id, x, y, z]

		if(strcmp(buffer, "end") == 0) {
			sender.sendMsg("ok", 2);
			break;
		}
		else { // it's a position
			if(parseList(buffer, pos)) {	
				if((int)pos[0] > max) {
					max = (int)pos[0];
				}

				neurons_.push_back(Neuron((int)pos[0], pos[1], pos[2], pos[3]));
			}
			else {
				if((int)pos[0] > max) {
					max = (int)pos[0];
				}

				neurons_.push_back(Neuron((int)pos[0], generateRandomNumber(RANDOM_POS_LOW, RANDOM_POS_HIGH), generateRandomNumber(RANDOM_POS_LOW, RANDOM_POS_HIGH), generateRandomNumber(RANDOM_POS_LOW, RANDOM_POS_HIGH)));
			}
			sender.sendMsg("ok", 2);
		}
	}

	// instead of using the vector, we dinamically create
	// an array of Neurons of size equal to the highest
	// id, so that we put every neuron in the array
	// according to its id. It will be simpler and
	// faster to acces them after.
	neurons = new Neuron[max + 1];

	for(int i = 0; i < neurons_.size(); i++) {
		neurons[neurons_.at(i).getId()] = neurons_.at(i);
	}
	nb_neurons = max + 1;
	// we don't waste space since we clear
	neurons_.clear();
}

void GraphicsSimulator::receive_connections() {
	char bufferParams[50];
	double params[3];
	
	char* bufferConn;
	int lengthConn;
	double* conn;
	int nbConn;
	
	while(true){
		listener.receiveMsg(bufferParams, 50); // [neuron id, nb connections, length of next msg]
		
		if(strcmp(bufferParams, "end") == 0) {
			sender.sendMsg("ok", 2);
			break;
		}
		else { // it's a connection
			if(parseList(bufferParams, params)) {
				nbConn = (int)params[1];
				lengthConn = (int)params[2];
				bufferConn = new char[lengthConn + 1];
				conn = new double[nbConn];
				
				sender.sendMsg("param_ok", 8);
				
				listener.receiveMsg(bufferConn, lengthConn + 1); // [target1, target2, ...]

				if(parseList(bufferConn, conn)) {
					for(int i = 0; i < nbConn; i++) {
						neurons[(int)params[0]].addConnection((int)conn[i]);
					}
					sender.sendMsg("msg_ok", 6);
				}
				else {
					printf("Problem receiving connections. Program will exit.");
					exit(1);
				}
				
				delete[] bufferConn;
				delete[] conn;
			}
			else {
				printf("Problem receiving connections. Program will exit.");
				exit(1);
			}
			sender.sendMsg("ok", 2);
		}
	}
}

// end of connection management












// simulation management
void GraphicsSimulator::start() {
	pthread_t spike_detector;
	
	// Variables for handling pause: basically, the current_time will just substract the total pause time
	long pausedTime = 0;
	long start_paused_time = 0;
	int eventType = EVENT_NOTHING;
	stopped = false;
	
	pthread_create(&spike_detector, NULL, detect_spikes, (void*)this);
	sender.sendMsg("simulate", 8);
	init_time = SDL_GetTicks();

	do
	{
		eventType = handleEvents();
		
		switch(eventType)
		{
		case EVENT_QUIT:
			sender.sendMsg("quit", 4);
			break;
		case EVENT_STEP_CHANGED:
			// We calculate the initial time so that the simulation would have the current percentage
			// with the new simulation step. This ensures that when the user changes the simulation speed,
			// the simulated time does not change.
			init_time = (SDL_GetTicks() - pausedTime) - (curr_time * simulation_step);
			break;
		case EVENT_STOP:
			start_paused_time = SDL_GetTicks();
			sender.sendMsg("stop", 4);
			stopped = true;
			break;
		case EVENT_RESUME:
			pausedTime += SDL_GetTicks() - start_paused_time;
			sender.sendMsg("resume", 6);
			stopped = false;
			break;
		}
		
		if (!stopped)
		{
			// The simulation_step variable divides the elapsed time and provides a way to slow down 
			// or accelerate the simulation (no more than real time)
		    curr_time = ((SDL_GetTicks() - pausedTime) - init_time) / simulation_step;		
		
		    if(curr_time <= sim_time) {
		    	for(int i=0; i < nb_neurons; i++) {
					neurons[i].update(curr_time);
				}
		    }
		}
		
		window.update(curr_time);
		camera.update();
		window.draw();
		
		SDL_Delay(SIMULATION_DELTA);
	} while(eventType != EVENT_QUIT);
	pthread_cancel(spike_detector);
}




void* detect_spikes(void* simulator_) {
	GraphicsSimulator* simulator = (GraphicsSimulator*) simulator_;
	char bufferParams[50];
	double params[2];
	
	char* bufferSpikes;
	int lengthSpikes;
	double* spikes;
	int nbSpikes;
	double spike_time;
	
	while(true){
		simulator->listener.receiveMsg(bufferParams, 50); // [nb spikes, length of next message]
		
		if(strcmp(bufferParams, "finish") == 0) {
			simulator->sender.sendMsg("ok", 2);
			break;
		}
		else { // it's a spike train
			if(parseList(bufferParams, params)) {
				nbSpikes = (int)params[0];
				lengthSpikes = (int)params[1];
				bufferSpikes = new char[lengthSpikes + 1];
				spikes = new double[nbSpikes];
				
				simulator->sender.sendMsg("ok", 2);
				
				simulator->listener.receiveMsg(bufferSpikes, lengthSpikes + 1); // [id1, id2, ..., time1, time2, ...]

				if(parseList(bufferSpikes, spikes)) {
					for(int i = 0; i < nbSpikes / 2; i++) {
						spike_time = spikes[i + nbSpikes / 2];

						// Checking this ensures that the spike list will be sorted in time
						if(spike_time >= simulator->curr_time) {
							simulator->neurons[(int)spikes[i]].fire(spike_time);
						}
					}
				}
				
				delete[] bufferSpikes;
				delete[] spikes;
			}
		}
	}
}





// events management


// Handles events and returns
// the event type
int GraphicsSimulator::handleEvents() {
	SDL_PollEvent(&event);
	
	switch(event.type)
	{
	case SDL_QUIT:
		return EVENT_QUIT;
		break;
	case SDL_KEYDOWN:
		switch(event.key.keysym.sym)
		{
		case SDLK_UP:
			camera.up();
			break;
		case SDLK_DOWN:
			camera.down();
			break;
		case SDLK_LEFT:
			camera.left();
			break;
		case SDLK_RIGHT:
			camera.right();
			break;
		case SDLK_KP_MINUS:
			if(!minus_pressed) {
				minus_pressed = true;
				
				if(!stopped) {
					incrementSimulationStep();
					return EVENT_STEP_CHANGED;
				}
			}
			break;
		case SDLK_KP_PLUS:
			if(!plus_pressed) {
				plus_pressed = true;
				
				if(!stopped) {
					decrementSimulationStep();
					return EVENT_STEP_CHANGED;
				}
			}
			break;	
		case SDLK_p:
			if(!p_pressed) {
				p_pressed = true;
				stopped = !stopped;
				if(!stopped) {
					return EVENT_RESUME;
				} else {
					return EVENT_STOP;
				}
			}
			break;
		case SDLK_w:
			camera.forward();
			break;
		case SDLK_s:
			camera.backward();
			break;
		case SDLK_1:
			if(!k1_pressed) {
				k1_pressed = true;
				camera.setMode(MODE_CENTERED);
			}
			break;
		case SDLK_2:
			if(!k2_pressed) {
				k2_pressed = true;
				camera.setMode(MODE_FREE);
			}
			break;
		}
		break;
	case SDL_KEYUP: // Events for making each button pressing non continuous
		switch(event.key.keysym.sym)
		{
		case SDLK_KP_MINUS:
			if(minus_pressed) {
				minus_pressed = false;
			}
			break;
		case SDLK_KP_PLUS:
			if(plus_pressed) {
				plus_pressed = false;
			}
			break;
		case SDLK_p:
			if(p_pressed) {
				p_pressed = false;
			}
			break;
		case SDLK_1:
			if(k1_pressed) {
				k1_pressed = false;
			}
			break;
		case SDLK_2:
			if(k2_pressed) {
				k2_pressed = false;
			}
			break;
		}
		break;
	case SDL_MOUSEBUTTONDOWN:
		switch(event.button.button) 
		{
		case SDL_BUTTON_WHEELUP:
			camera.forward();
			break;
		case SDL_BUTTON_WHEELDOWN:
			camera.backward();
			break;
		}
		break;
	}
	
	return EVENT_NOTHING;
}



// We decrement or increment in a way that allows
// the user to do it faster. It actually calculates
// the order of the current simulation_step
// and adds or subtracts a value of this same order.

void GraphicsSimulator::incrementSimulationStep() {
	int tmp = simulation_step;
	int order = 0;
	
	while(tmp >= 1) {
		tmp /= 10;
		order++;
	}
	int inc = pow(10, order - 1);
	
	simulation_step += inc;
	if(simulation_step > HIGH_BOUND_SIM_STEP ) {
		simulation_step = HIGH_BOUND_SIM_STEP;
	}
}



void GraphicsSimulator::decrementSimulationStep() {
	int tmp = simulation_step;
	int order = 0;
	
	while(tmp >= 1) {
		tmp /= 10;
		order++;
	}
	int inc = pow(10, order - 1);
	
	if(simulation_step - inc > 0) {
		simulation_step -= inc;
	} else {
		simulation_step -= inc / 10;
	}
	
	if(simulation_step < LOW_BOUND_SIM_STEP ) {
		simulation_step = LOW_BOUND_SIM_STEP;
	}
}




