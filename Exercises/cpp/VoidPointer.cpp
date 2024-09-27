#include <bits/stdc++.h>
using namespace std;

void funct(void * arg, int N){
	char* c = (char*) arg;
	for(int i = 0; i < N; i++){
		cout << *c;
		c++;
	}
	cout << '\n';
}

int main(){
	char str[] = "Hello, world!";
	funct(str, 14);
}
