#ifndef INCLUDED_HTTPSERVER
#define INCLUDED_HTTPSERVER

#include<arpa/inet.h>
#include<string>
#include <sqlite3.h>

namespace http {
    class TcpServer {
		public:
			TcpServer(std::string ip_address, int port);
			~TcpServer();
			void startListen();

		private:
			std::string m_ip_address;
			int m_port;
			int m_socket;
			int m_new_socket;
			long m_incoming_message;
			sockaddr_in m_socketAddress;
			unsigned int m_socketAddress_len;
			std::string m_serverMessage;

			sqlite3* db;

			int startServer();
			void closeServer();
			std::string buildResponse(std::string);
			void acceptConnection(int &new_socket);
			void sendResponse();
    };
}

#endif
