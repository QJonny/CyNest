#!/bin/sh

cwd= pwd

cd $1
g++ main.cpp graphics_simulator.cpp window.cpp network.cpp socket.cpp tools.cpp `sdl-config --cflags --libs` -lGLU -lSDL -lpthread -lSDL_ttf -o g_simulator
	
cd $cwd
