#include <iostream>

#include "src/LBFGS.h"


float fn(float* in) {
    return in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
}

void fn_grad(float* in, float* out) {
    out[0] = 2*in[0];
    out[1] = 2*in[1];
    out[2] = 2*in[2];
}

int main() {
    float in[3];
    in[0] = 1.0;
    in[1] = 1.0;
    in[2] = 2.0;
    printf("f(%f, %f, %f) = %f\n", in[0], in[1], in[2], fn(in));
    float out[3];
    fn_grad(in, out);
    printf("grad of f(%f, %f, %f) = [ %f, %f, %f ]\n", in[0], in[1], in[2], out[0], out[1], out[2]);



    return 0;
}
