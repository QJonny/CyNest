#include "socket.h"

using namespace std;

void Socket::init(int port_, char* host_) {
	#if _WIN32
		WSADATA wsaData;
		WORD version;
		version = MAKEWORD( 2, 0 );
		WSAStartup( version, &wsaData );
		WSACleanup();
		recsize = sizeof(ssin);
		crecsize = sizeof(csin);
	#else
		recsize = sizeof(ssin);
		crecsize = sizeof(csin);
	#endif
	port = port_;
	sock = socket(AF_INET, SOCK_STREAM, 0);
	ssin.sin_addr.s_addr = htonl(INADDR_ANY);
	ssin.sin_family = AF_INET;
	ssin.sin_port = htons(port_);
}
	
	
void Socket::acceptConnection(int port_, char* host_) {
	init(port_, host_);
	bind(sock, (sockaddr*)&ssin, recsize);
	listen(sock, 1);
	csock = accept(sock, (sockaddr*)&csin, &crecsize);
	listener = 1;
}

void Socket::initiateConnection(int port_, char* host_) {
	init(port_, host_);
	connect(sock, (sockaddr*)&ssin, recsize);
	listener = 0;
}


void Socket::sendMsg(char* msg, int length) {
	send(sock, msg, length, 0);
}

void Socket::receiveMsg(char* buffer, int length) {
	memset(buffer, 0, length);
	recv(csock, buffer, length, 0);
}

void Socket::destroy(){
	shutdown(sock, 2);
	if(listener = 1){
		shutdown(csock, 2);
	}
}
	
	
