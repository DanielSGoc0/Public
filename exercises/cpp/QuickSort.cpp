/*
This is an exercise in creating references
and dereferencing. It's a quick sort algorithm.
*/
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <list>
using namespace std;

bool x_axis = false;

class Point {
	private:
		int x, y;

	public:
		Point(int a, int b){
			x = a;
			y = b;
		}

		void transpose();

		int get_x() const {
			return x;
		}

		int get_y() const {
			return y;
		}
};

void Point::transpose(){
	swap(x, y);
}

ostream& operator<<(ostream& out, Point& point){
	out << "(" << point.get_x() << ", " << point.get_y() << ")";
	return out;
}

template<typename T, size_t N> void print_array(ostream& out, T (&arr)[N]){
	for(auto i = 0; i < N; i++){
		out << arr[i] << ' ';
	}
	out << '\n';
}

template<typename T> void print_interval(ostream& out, typename T::iterator b, typename T::iterator e){
	for(; b != e; b++){
		out << *b << ' ';
	}
	out << '\n';
}

bool operator<(const Point& P, const Point& Q){
	if(x_axis){
		if(P.get_x() < Q.get_x()){
			return true;
		} else if(Q.get_x() < P.get_x()){
			return false;
		} else {
			return P.get_y() < Q.get_y();
		}
	} else {
		if(P.get_y() < Q.get_y()){
			return true;
		} else if(Q.get_y() < P.get_y()){
			return false;
		} else {
			return P.get_x() < Q.get_x();
		}
	}
}


// QuickSort based on templates and bidirectional iterators.
template<typename T> void QuickSort(typename T::iterator b, typename T::iterator e){
	if(b == e || next(b) == e){
		return;
	}

	// pivot is *b
	typename T::iterator i = next(b);
	typename T::iterator j = prev(e);

	while(i != j){
		if(*i < *b){
			i = next(i);
			continue;
		} else if(!(*j < *b)){
			j = prev(j);
			continue;
		}
		swap(*i, *j);
		j = prev(j);
	}
	if(*i < *b) swap(*b, *i);

	QuickSort<T>(b, i);
	QuickSort<T>(i, e);
}


// QuickSort based on a reference to array and interval. The comparison operator < must overriden.
template<typename T> void QuickSort(T arr[], int b, int e){
	if(b == e || b + 1 == e){
		return;
	}

	// pivot is *b
	int i = b + 1;
	int j = e - 1;

	while(i != j){
		if(arr[i] < arr[b]){
			i++;
			continue;
		} else if(!(arr[j] < arr[b])){
			j--;
			continue;
		}
		swap(arr[i], arr[j]);
		j--;
	}
	if(arr[i] < arr[b]) swap(arr[b], arr[i]);

	QuickSort<T>(arr, b, i);
	QuickSort<T>(arr, i, e);
}


// append N random points at the end of vector.
void random_points(vector<Point>& points, int N){
	for(int i = 0; i < N; i++){
		points.push_back(Point(rand(), rand()));
	}
}

// save results to a file.
void save_to_output_file(string fname, vector<Point>& points){
	ofstream out;
	out.open("./" + fname);
	for(auto it = begin(points); it != end(points); it++){
		out << *it << ' ';
	}
	out.close();
}

// append points from file at the end of vector.
void get_input(string fname, vector<Point>& points){
	ifstream in;
	in.open("./" + fname, ios::in);
	in >> ws;
	char c;
	while(!in.eof() && points.size() < 12){
		int x, y;
		in.get(c);
		in >> x;
		in >> ws;
		in.get(c);
		in >> ws;
		in >> y;
		in >> ws;
		in.get(c);
		in >> ws;
		points.push_back(Point(x, y));
	}
	in.close();
}

int main(){
	vector<Point> points;
	get_input("in", points);

	x_axis = true;
	QuickSort<vector<Point> >(begin(points), end(points));
	save_to_output_file("out", points);

	x_axis = false;
	list<Point> plist;
	for(auto it = begin(points); it != end(points); it++){
		plist.push_front(*it);
	}
	QuickSort<list<Point> >(begin(plist), end(plist));
	for(auto it = begin(plist); it != end(plist); it++){
		cout << *it << endl;
	}
}
