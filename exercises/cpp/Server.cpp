/*
idea for a project:
user sends an HTTP requests (via postman) to add/change/extract
information from SQLite database. All communication is handled by
this c++ program. This is an exercise in communication.
*/

/*
There are 3 types of requests:
GET: returns all contents of database.
POST (val): append the value val at the end of database.
DELETE (val): delete all the values val from the database.
*/

// Much of the code was copied from:
// https://osasazamegbe.medium.com/showing-building-an-http-server-from-scratch-in-c-2da7c0db6cb7


#include "HTTPServer.h"
#include <iostream>
using namespace std;
using namespace http;


// int callback(void *NotUsed, int argc, char **argv, char **azColName){
// 	for(int i = 0; i < argc; i++){
// 		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
// 	}

// 	printf("\n");
// 	return 0;
// }

int main(){
	
	TcpServer server = TcpServer("0.0.0.0", 8080);
	server.startListen();

	return 0;
}

// int main(){

// 	sqlite3* db;
// 	char* ErrMsg = 0;

// 	sqlite3_open("test.db", &db);

// 	// sqlite3_exec(db, "SELECT * FROM sqlite_master;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;
// 	// sqlite3_exec(db, "CREATE TABLE vals(val INTEGER);", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;

// 	// sqlite3_exec(db, "INSERT INTO vals(val) VALUES (3), (3), (7), (8);", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;

// 	// sqlite3_exec(db, "DELETE FROM vals WHERE val < 9 LIMIT 1;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;

// 	// sqlite3_exec(db, "SELECT val FROM vals;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;

// 	// sqlite3_exec(db, "SELECT * FROM sqlite_master;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;
// 	// sqlite3_exec(db, "DROP TABLE tbl1;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;
// 	// sqlite3_exec(db, "SELECT * FROM sqlite_master;", callback, 0, &ErrMsg);
// 	// cout << ErrMsg << endl;

	
// 	//  stmt = sqlite3_prepare16_v3(db, ".tables", ErrMsg)

// 	cout << sqlite3_close(db) << endl;
// }
