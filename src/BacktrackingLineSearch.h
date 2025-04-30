#ifndef BACKTRACKINGLINESEARCH_H
#define BACKTRACKINGLINESEARCH_H

#include <functional>

using namespace std;

template <typename T>

class BacktrackingLineSearch {
public:
	BacktrackingLineSearch(T* x, T* g, T* q, int dof, function<T(T*)> user_func) : X(x), G(g), DOF(dof), func(user_func) {
		x_plus_step = (T*)malloc(DOF * sizeof(T));
		search_direction = (T*)malloc(DOF * sizeof(T));
		for (int i = 0; i < DOF; i++) {
			search_direction[i] = -q[i];
		}
	
	}
	T linesearch() {
		int max_it = 1000;
		T c = 0.5;
		T tau = 0.75;
		T m = dot_product(G, search_direction, DOF);
		T step_size = 1; 
		for (int i = 0; i < max_it; i++) {
			for (int j = 0; j < DOF; ++j) {
				x_plus_step[j] = X[j] + (step_size * search_direction[j]);
			}
			if (func(X) - func(x_plus_step) >= step_size * -c * m) {
				return step_size;
			}
			else {
				step_size *= tau;
			}
		}
		return 0;
	}
	~BacktrackingLineSearch() {
		free(x_plus_step);
		free(search_direction);
	}
private:
	T* X;
	T* G;
	T* x_plus_step;
	T* search_direction;
	int DOF;
	function<T(T*)> func;

	T dot_product(T* a, T* b, int n) {
		T result = 0;
		for (int i = 0; i < n; ++i) {
			result += a[i] * b[i];
		}
		return result;
	}

};

#endif