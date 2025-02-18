//
// Created by Matthew Speranza on 2/13/25.
//

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <cstddef>
#include <type_traits>

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
    virtual T* minimize_step(T* X, T* G, T step_size) = 0;
};

#endif //MINIMIZER_H
