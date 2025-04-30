#include <iostream>
#include <cstdlib>


template <typename T>
class SteepestDescent : Minimizer<T> {
public:
	SteepestDescent(int DOF) : DOF(DOF) {
		prev_positions = (T*)malloc(DOF * sizeof(T));
		prev_gradient = (T*)malloc(DOF * sizeof(T));
		max_step_size = 1.0;
	}

	void init(T* X, T* G, T steepest_descent_step_size) override {
		for (int i = 0; i < DOF; ++i) {
			prev_positions[i] = X[i];
			prev_gradient[i] = G[i];
		}
		for (int i = 0; i < DOF; ++i) {
			X[i] -= (steepest_descent_step_size * G[i]);
		}
	}
	void minimize() override {
	}

	void minimize_step(T* X, T* G) override {
		T step_size = Minimizer::adaptive_step_size(X, prev_positions, G, prev_gradient, DOF);
		if (step_size > max_step_size) {
			step_size = 0.1;
		}
		std::cout << "step size: " << step_size << "\n";
		for (int i = 0; i < DOF; ++i) {
			prev_positions[i] = X[i];
			prev_gradient[i] = G[i];
		}
		for (int i = 0; i < DOF; ++i) {
			X[i] -= (step_size * G[i]);
		}
	}

private:
	int DOF;
	T max_step_size;
	T* prev_positions;
	T* prev_gradient;
};