
#ifndef LBFGS_H
#define LBFGS_H
#include <type_traits>

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <type_traits>

#include "Minimizer.h"

template <typename T>
class LBFGS : public Minimizer<T> {
    static_assert(std::is_floating_point<T>::value, "L-BFGS can only be used with floating point types");

public:
    LBFGS(int m, int DOF) : m(m), DOF(DOF) {
        min_step_size = sizeof(T) == 4 ? 1e-4 : 1e-7;
        rho = malloc(m*sizeof(T));
        alpha = malloc(m*sizeof(T));

        search = malloc(DOF*sizeof(T));
        prev_positions = malloc(DOF*sizeof(T));
        prev_gradient = malloc(DOF*sizeof(T));

        s = malloc(m*DOF*sizeof(T));
        y = malloc(m*DOF*sizeof(T));
    }

    /**
     * TODO: Free all memory used by LBFGS
     */
    ~LBFGS() override {
    }

    /**
     * TODO: Implement steepest decent step
     * Sets variables and returns 1
     * @param X Initial Position
     * @param G Initial Position Grad
     */
    void init(T* X, T* G, T steepest_decent_step_size) override {
    }

    /**
     * This will loop over minimize step and a line search. Implement later
     */
    void minimize() override {
    }

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
    T* minimize_step(T* X, T* G, T step_size) override {

        return X;
    }

    /**
     * TODO: Update variables s, y, rho, alpha, prev_positions, prev_gradients
     */
    void update() {

    }

private:
    int m; // Number of previous gradients to use for hessian approximation (5-7)
    int DOF; // Degrees of freedom
    T min_step_size; // terminate minimization with steps smaller than this number

    T gamma; // s projected onto y
    T *rho; // [m] rho_{i} = (s_{i}^T * y_{i}
    T *alpha; // [m] alpha_{i} = rho_{i} * s_{i}^T * y_{i}

    T *X;
    T *G;
    T *search; // [DOF] L-BFGS search direction
    T *prev_positions; // [DOF] x_{i-1}
    T *prev_gradient; // [DOF] g_{i-1}

    T *s; // [m*DOF] x_{i+1} - x_{i} = s_{i}
    T *y; // [m*DOF] grad_{i+1} - grad_{i} = y_{i}
};

#endif //LBFGS_H
