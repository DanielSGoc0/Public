/*
A custom class to be used in python exercise.
*/
#include <iterator>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace std;
namespace py = pybind11;


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


class MyClass{
	public:
		vector<double> MyVec;

		MyClass(vector<double> MyVec_in){
			for(auto it = begin(MyVec_in); it != end(MyVec_in); it++){
				MyVec.push_back(*it);
			}
		}

		void run(){
			QuickSort<vector<double> >(begin(MyVec), end(MyVec));
		}

		vector<double> getVec(){
			return MyVec;
		}
};

PYBIND11_MODULE(QuickSort, m) {
    py::class_<MyClass>(m, "MyClass")
        .def(py::init<const std::vector<double> >(), "initialize with a vector of real numbers")
        .def("run", &MyClass::run, "this function sorts your array.")
        .def("getVec", &MyClass::getVec, "returns the array.");
}
