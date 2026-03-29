#include <bits/stdc++.h>
using namespace std;

const int N = 1024;
const int INTERVAL = 10;
const int LEVELS = 13;

double sqrts[N];
double coeffs[N];
vector<double> xs;

default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);


void initialize(){
	generator.seed(time(0));

	for(int n = 0; n < N; n++){
		sqrts[n] = 1.0/sqrt(n + 1);
	}

	for(int l = 0; l < LEVELS; l++){
		for(int j = 0; j < (1<<l); j++){
			double x = INTERVAL * (-1 + (2*j + 1.0) / (1<<l));
			xs.push_back(x);
		}
	}
}

void new_random_function(){
	for(int n = 0; n < N; n++){
		coeffs[n] = distribution(generator);
	}

	if(coeffs[0] < 0){
		coeffs[0] = -coeffs[0];	
	}
}

double f(double x){
	double y = 0;
	for(int n = N - 1; n >= 0; n--){
		y *= x * sqrts[n];
		y += coeffs[n];
	}
	return y * exp(-x*x/2.0);
}

vector<double> new_point(){
	while(true){
		new_random_function();

		bool POSITIVE = true;

		for(int i = 0; i < xs.size(); i++){
			double y = f(xs[i]);
			if(y < 0){
				POSITIVE = false;
				break;
			}
		}
		if(POSITIVE){
			vector<double> RESULT;
			RESULT.push_back(coeffs[0]);
			RESULT.push_back(coeffs[1]);
			RESULT.push_back(2*coeffs[2] - coeffs[0]);
			return RESULT;
		}
	}
}

void add_new_points(int MAX){
	ofstream output_file;
	output_file.open("samples/N" + to_string(N) + "INTERVAL" + to_string(INTERVAL) + "LEVELS" + to_string(LEVELS) + ".txt");

	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::showpoint);
	cout.precision(12);
	output_file.setf(ios::fixed, ios::floatfield);
	output_file.setf(ios::showpoint);
	output_file.precision(12);

	for(int i = 0; i < MAX; i++){
		vector<double> res = new_point();

		output_file << res[0] << " " << res[1] << " " << res[2] << "\n";
		cout << i << "   \t" << res[0] << "   \t" << res[1] << "   \t" << res[2] << "\n";
	}
}

int main(){
	initialize();

	add_new_points(10000);
}
