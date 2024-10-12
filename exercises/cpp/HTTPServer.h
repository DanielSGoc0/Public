#include "HTTPServerClass.h"

#include <sys/socket.h>
#include <stdexcept>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sqlite3.h>
#include <cstring>

namespace http {
	std::string response;

	int callbackGET(void *NotUsed, int argc, char **argv, char **azColName){
		for(int i = 0; i < argc; i++){
			response.append(argv[i]);
			response.append(", ");
		}

		return 0;
	}

	int callbackPOST(void *NotUsed, int argc, char **argv, char **azColName){
		return 0;
	}

	int callbackDELETE(void *NotUsed, int argc, char **argv, char **azColName){
		return 0;
	}


	const int BUFFER_SIZE = 30720;

	TcpServer::TcpServer(std::string ip_address, int port):
		m_ip_address(ip_address), m_port(port), m_socket(), m_new_socket(),
		m_incoming_message(), m_socketAddress(), m_socketAddress_len(sizeof(m_socketAddress)),
		m_serverMessage(buildResponse("Hello from your Server :)")) {

		std::cout << "Constructor called" << std::endl;

		m_socketAddress.sin_family = AF_INET;
		m_socketAddress.sin_port = htons(m_port);
		m_socketAddress.sin_addr.s_addr = inet_addr(m_ip_address.c_str());
		startServer();
	};

	TcpServer::~TcpServer(){
		
		std::cout << "Destructor called" << std::endl;

		closeServer();
	};

	int TcpServer::startServer() {
		m_socket = socket(AF_INET, SOCK_STREAM, 0);
		if (m_socket < 0) {
			throw std::invalid_argument("Cannot create socket");
			return 1;
		}
		if (bind(m_socket,(sockaddr *)&m_socketAddress, m_socketAddress_len) < 0) {
			throw std::invalid_argument("Cannot connect socket to address");
			return 1;
		}

		sqlite3_open("test.db", &db);

		std::cout << "successfully started server..." << std::endl;

		return 0;
	}

	void TcpServer::closeServer() {
		close(m_socket);
		close(m_new_socket);
		sqlite3_close(db);
		exit(0);
	}

	void TcpServer::startListen() {
		if (listen(m_socket, 20) < 0) {
			throw std::invalid_argument("Socket listen failed");
		}

		std::cout << "\n*** Listening on ADDRESS: " 
			<< inet_ntoa(m_socketAddress.sin_addr) 
			<< " PORT: " << ntohs(m_socketAddress.sin_port) 
			<< " ***\n\n";

        int bytesReceived;

        while (true) {
            std::cout << "====== Waiting for a new connection ======\n\n\n";
            acceptConnection(m_new_socket);

            char buffer[BUFFER_SIZE] = {0};
            bytesReceived = read(m_new_socket, buffer, BUFFER_SIZE);
            if (bytesReceived < 0) {
                std::cout << "Failed to read bytes from client socket connection" << std::endl;
				exit(1);
            }

			std::cout << "------ Received Request from client ------\n\n";

			std::istringstream in_stream(buffer);
			std::string request_type;
			in_stream >> request_type;
			int val = -1;
			char* ErrMsg;

			char c;
			for(int i = 0; i < 8; i++){
				in_stream >> c;
			}
			in_stream >> val;

			std::cout << "querying database... " << request_type << "  " << val << std::endl;

			if(request_type[0] == 'G'){
				response.clear();
				sqlite3_exec(db, "SELECT val FROM vals;", callbackGET, 0, &ErrMsg);
			} else if(request_type[0] == 'P'){
				response = "added new value!";
				sqlite3_exec(db, ("INSERT INTO vals(val) VALUES (" + std::to_string(val) + ");").c_str(), callbackPOST, 0, &ErrMsg);
			} else if(request_type[0] == 'D'){
				response = "removed values!";
				sqlite3_exec(db, ("DELETE FROM vals WHERE val = " + std::to_string(val) + ";").c_str(), callbackDELETE, 0, &ErrMsg);
			} else {
				response = "invalid request";
			}

			m_serverMessage = buildResponse(response);

			// std::cout << ErrMsg << std::endl;

			std::cout << "sending response..." << std::endl;

            sendResponse();
            close(m_new_socket);
        }
	}

	void TcpServer::acceptConnection(int &new_socket) {
		new_socket = accept(m_socket, (sockaddr *)&m_socketAddress, &m_socketAddress_len);

		if (new_socket < 0) {
			std::cout << "Server failed to accept incoming connection from ADDRESS: " 
				<< inet_ntoa(m_socketAddress.sin_addr) << "; PORT: " 
				<< ntohs(m_socketAddress.sin_port) << std::endl;
		}
	}

	std::string TcpServer::buildResponse(std::string content) {
        std::string htmlFile = "<!DOCTYPE html><html lang=\"en\"><body><h1> HOME </h1><p> " + content + " </p></body></html>";
        std::ostringstream ss;
        ss << "HTTP/1.1 200 OK\nContent-Type: text/html\nContent-Length: " << htmlFile.size() << "\n\n"
           << htmlFile;

        return ss.str();
    }

	void TcpServer::sendResponse() {
        long bytesSent;

		// std::cout << m_serverMessage << std::endl;

        bytesSent = write(m_new_socket, m_serverMessage.c_str(), m_serverMessage.size());

        if (bytesSent == m_serverMessage.size()) {
            std::cout << "------ Server Response sent to client ------\n\n";
        } else {
            std::cout << "Error sending response to client" << std::endl;
        }
    }
}
