
#ifndef LBFGS_H
#define LBFGS_H
#include <type_traits>

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <functional>
#include "BacktrackingLineSearch.h"
#include "Minimizer.h"

template <typename T>
class LBFGS : public Minimizer<T> {
    static_assert(std::is_floating_point<T>::value, "L-BFGS can only be used with floating point types");

public:
    LBFGS(int m, int DOF, std::function<T(T*)> user_func, std::function<void(T*, T*)> user_grad) : m(m), DOF(DOF), func(user_func), grad(user_grad) {
        minimized = false;
        min_step_size = sizeof(T) == 4 ? 1e-4 : 1e-7;
        max_step_size = 1.0;
        rho = (T*) calloc(m, sizeof(T));
        alpha = (T*) malloc(m*sizeof(T));
        gamma = (T*) malloc(sizeof(T));
        beta = (T*) malloc(sizeof(T));
        step_size = 0.0;
            

        G = (T*) malloc(DOF*sizeof(T));
        X = (T*) malloc(DOF*sizeof(T));

        search = (T*) malloc(DOF*sizeof(T));
        prev_positions = (T*)malloc(DOF * sizeof(T));
        prev_gradient = (T*)malloc(DOF * sizeof(T));

        q = (T*) malloc(DOF*sizeof(T));
        s = (T*) calloc(m*DOF, sizeof(T));
        y = (T*) calloc(m*DOF, sizeof(T));
    }

    /**
     * TODO: Free all memory used by LBFGS
     */
    ~LBFGS() override {
        free(rho);
        free(alpha);
        free(gamma);
        free(beta);
        free(G);
        free(X);
        free(search);
        free(prev_positions);
        free(prev_gradient);
        free(q);
        free(s);
        free(y);
    }

/**
 * TODO: Implement steepest decent step
 * Sets variables and returns 1
 * @param X Initial Position
 * @param G Initial Position Grad
 */
void init(T* X, T* G, T steepest_descent_step_size) override {
    for (int i = 0; i < DOF; ++i) {
        prev_positions[i] = X[i];
        prev_gradient[i] = G[i];
    }

    for (int i = 0; i < DOF; ++i) {
        X[i] -= steepest_descent_step_size * G[i];
        s[i + (m - 1) * DOF] = X[i] - prev_positions[i];
        y[i + (m - 1) * DOF] = G[i] - prev_gradient[i];
        this->X[i] = X[i];
    }
    grad(this->X, this->G);
}

/**
 * This will loop over minimize step and a line search. Implement later
 */
void minimize() override {
    bool found_minimum = false;
    int j = 0;
    int maxcount2 = 1000000;
    while (!found_minimum && j < maxcount2) {
        ++j;
        minimize_step(X, G);
        printf("f(");
        for (int i = 0; i < DOF; i++) {
            if (i != DOF - 1) {
                printf("%f, ", X[i]);
            }
            else {
                printf("%f) = %f\n", X[i], func(X));
            }
        }
        grad(X, G);
        printf("grad of f(");
        for (int i = 0; i < DOF; i++) {
            if (i != DOF - 1) {
                printf("%f, ", X[i]);
            }
            else {
                printf("%f) = [", X[i]);
            }
        }
        for (int i = 0; i < DOF; i++) {
            if (i != DOF - 1) {
                printf("%f, ", G[i]);
            }
            else {
                printf("%f]\n", G[i]);
            }
        }
        double sum_of_residuals = 0;
        for (int i = 0; i < DOF; ++i) {
            sum_of_residuals += G[i] * G[i];

        }
        sum_of_residuals = sqrtf((sum_of_residuals / DOF));
        if (is_minimized() == true || sum_of_residuals < 1e-7) {
            std::cout << "found minimum\n";
            std::cout << "Number of Steps: " << j << "\n";
            found_minimum = true;
        }

    }
}
T line_search(T(*function)(T*), void(*gradient)(T*, T*)) {
    T step = 1;
    T* x_plus_step = (T*)malloc(DOF * sizeof(T));
    T c1 = 1e-3;

    for (int i = 0; i < 1000; ++i) {
        for (int i = 0; i < DOF; ++i) {
            x_plus_step[i] = X[i] + (step * q[i]);
        }
        if (function(x_plus_step) <= function(X) + (c1 * step * dot_product(q, gradient(X), DOF))) {
            free(x_plus_step);
            return step;
        }
        else {
            step *= 0.5;
        }
    }
    free(x_plus_step);
    return 0;
}


// checks if a large enough step is taken
/*
bool armijo_rule(T* function(T*), T* gradient(T*), step) {
    T *x_plus_step = (T*) malloc(DOF*sizeof(T));
    for (int i = 0; i < DOF; ++i) {
        x_plus_step[i] = X[i] + (step * q[i]);
    }
    bool result = (function(x_plus_step) <= function(X) + (c1 * step * dot_product(q, gradient(X), DOF);
    free(x_plus_step);
    return result;
}*/
/*bool curvature_condition(T* function(T*), T* gradient(T*), step) {
    T *x_plus_step = (T*)malloc(DOF * sizeof(T));
    for (int i = 0; i < DOF; ++i) {
        x_plus_step[i] = X[i] + (step * q[i]);
    }
    bool result = (abs(dot_product(q, gradient(x_plus_step), DOF)) <= c2 * abs(dot_product(q, gradient(X);// abs value makes it strong wolfe condtions
    free(x_plus_step);
    return result;
}*/
/**
 * TODO: Implement the below algorithm and call the update() method to store information gained from X and G
 * TODO: Return x_i-1 + step_size * search
 * 1. Compute new L-BFGS step direction
 *   Pseudocode from wikipedia:
 *   q = g.i // search direction to be updated
 *   for j = i-1 to i-m:
 *     alpha.j = rho.j * s.j.T * q // dot product
 *     q = q - alpha.j * y.j // vector scale & subtraction
 *   gamma.i = s.i-1.T * y.i-1 / y.i-1.T * y.i-1 // dot products in numerator and denominator
 *   q = gamma.i * q
 *   for j = i-m to i-1:
 *     beta = rho.j * y.j.T * q // dot product
 *     q = q + (alpha.j - beta) * s.j // vector scale & addition
 *   q = -q  // negate applied above instead of here most likely
 *   gamma = s.i.T * y.i / y.i.T * y.i
 *   rho.j = 1 / (y.j.T * s.j)
 */

void minimize_step(T* X_input, T* G_input) override {
    for (int i = 0; i < DOF; ++i) {
        X[i] = X_input[i];
        G[i] = G_input[i];
    }
    for (int i = 0; i < DOF; ++i) {
        q[i] = G[i];
    }
    //Minimizer::adaptive_step_size(X, prev_positions, G, prev_gradient, DOF);
    /*if (step_size > max_step_size) {
        step_size = 0.1;
    }*/
    //std::cout << "step size: " << step_size << "\n";
    update();
    for (int i = m - 1; i >= 0; --i) {
        alpha[i] = rho[i] * dot_product(s + i * DOF, q, DOF);


        for (int j = 0; j < DOF; ++j) {
            q[j] -= alpha[i] * y[i * DOF + j];
        }
    }

    T gamma = dot_product(s + (m - 1) * DOF, y + (m - 1) * DOF, DOF) / dot_product(y + (m - 1) * DOF, y + (m - 1) * DOF, DOF);
    for (int j = 0; j < DOF; ++j) {
        q[j] *= gamma;
    }

    for (int i = 0; i < m; ++i) {
        T beta = rho[i] * dot_product(y + i * DOF, q, DOF);

        for (int j = 0; j < DOF; ++j) {
            q[j] += (alpha[i] - beta) * s[i * DOF + j];
        }
    }
    BacktrackingLineSearch<double> backtracking(X, G, q, DOF, func);
    step_size = backtracking.linesearch(); 
    /*if (step_size <= 1e-7) {
        for (int i = 0; i < m * DOF; i++) {
            s[i] = 0;
            y[i] = 0;
        }
        for (int i = 0; i < m; i++) {
            rho[i] = 0;
            alpha[i] = 0;
        }
        LBFGS::init(X, G, 0.0001);  
    }*/
    std::cout << "step size: " << step_size << "\n";
    for (int j = 0; j < DOF; ++j) {
        X[j] += (step_size * -q[j]);
        X_input[j] = X[j];
    }
    
}


void update() {
    for (int i = 0; i < ((m - 1) * DOF); ++i) {
        s[i] = s[i + DOF];
        y[i] = y[i + DOF];

    }
    for (int i = 0; i < DOF; ++i) {
        s[(m - 1) * DOF + i] = X[i] - prev_positions[i];
        y[(m - 1) * DOF + i] = G[i] - prev_gradient[i];
    }

    for (int i = 0; i < m - 1; ++i) {
        rho[i] = rho[i + 1];
    }
    double s_dot_y = dot_product((s + (m - 1) * DOF), (y + (m - 1) * DOF), DOF);
    if (s_dot_y == 0) {
        std::cout << "Error: Dividing by zero. Function is likely already at a minimum.\n";
        minimized = true;
    }
 
    else {

        rho[m - 1] = 1 / s_dot_y;
    }
 

    for (int i = 0; i < DOF; ++i) {
        prev_positions[i] = X[i];
        prev_gradient[i] = G[i];
    }
      
}

bool is_minimized() {
    return minimized;
}
T get_step_size() {
    return step_size;
}

private:
    int m; // Number of previous gradients to use for hessian approximation (5-7)
    int DOF; // Degrees of freedom
    T min_step_size; // terminate minimization with steps smaller than this number
    T max_step_size; // ensures that lbfgs doesn't overshoot minimum
    bool minimized;

    T *beta; 
    T *gamma; // s projected onto y
    T *rho; // [m] rho_{i} = (s_{i}^T * y_{i}
    T *alpha; // [m] alpha_{i} = rho_{i} * s_{i}^T * y_{i}

    T *X;
    T *G;
    T *search; // [DOF] L-BFGS search direction
    T *prev_positions; // [DOF] x_{i-1}
    T *prev_gradient; // [DOF] g_{i-1}

    T* q;
    T *s; // [m*DOF] x_{i+1} - x_{i} = s_{i}
    T *y; // [m*DOF] grad_{i+1} - grad_{i} = y_{i}
    T step_size;

    std::function<T(T*)> func;
    std::function<void(T*, T*)> grad;

    T dot_product(T* a, T* b, int n) {
        T result = 0;
        for (int i = 0; i < n; ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
};

#endif //LBFGS_H
