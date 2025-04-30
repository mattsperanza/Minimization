//
// Created by Matthew Speranza on 2/13/25.
//

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <cstddef>
#include <type_traits>
#include <cmath>

template <typename T>
class Minimizer {
    static_assert(std::is_floating_point<T>::value, "Minimizer can only be used with floating point types");

public:
    /**
     * TODO: Free all memory associated with Minimizer not already freed by subclasses
     */
    virtual ~Minimizer() = default;
    virtual void init(T* X, T* G, T step_size) = 0;
    virtual void minimize() = 0;
    virtual void minimize_step(T* X, T* G) = 0;
    T adaptive_step_size(T* curr_pos,T* prev_pos,T* curr_grad,T* prev_grad, int DOF) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (int i = 0; i < DOF; ++i) {
            numerator += abs((curr_pos[i] - prev_pos[i]) * (curr_grad[i] - prev_grad[i]));
            denominator += (curr_grad[i] - prev_grad[i]) * (curr_grad[i] - prev_grad[i]);
        }
        return numerator / denominator;
    }
};

#endif //MINIMIZER_H
