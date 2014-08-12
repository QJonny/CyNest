#ifndef SOCKET_H
#define SOCKET_H

#include "headers.h"


#if defined (_WIN32)
    #include <winsock.h>
#else
    #include <sys/socket.h>
    #include <netinet/in.h>
    #include <arpa/inet.h>
    #include <sys/types.h>
#endif




class Socket {
private:
	#if _WIN32
		int sock;
		int csock;
		int recsize;
		int crecsize;
	#else
		socklen_t sock;
		socklen_t csock;
		socklen_t recsize;
		socklen_t crecsize;
	#endif

		int port;
		sockaddr_in ssin;

		sockaddr_in csin;
		int listener;


	void init(int port_, char* host_);
	
public:
	void acceptConnection(int port_, char* host_);
	
	void initiateConnection(int port_, char* host_);


	void sendMsg(char* msg, int length);
	
	void receiveMsg(char* buffer, int length);

	void destroy();
};


#endif
