#include <iostream>
#include <cmath>
#include "src/LBFGS.h"
#include "src/SteepestDescent.h"
#include "src/BacktrackingLineSearch.h"


double fn(double* in) {
    return (in[0] * in[0]) + (in[1] * in[1]);//(-2 * exp(-((in[0] + 2) * (in[0] + 2))/2)) - exp(-((in[0] - 2) * (in[0] - 2)) / 2) + ((in[0] * in[0]) / 20) + ((in[0] * in[1])/10) + exp(in[1] * in[1] * in[1] * in[1] * in [1]) + exp(-in[1] * in[1] * in[1] * in[1] * in[1]);
}

void fn_grad(double* in, double* out) {
    out[0] = (2 * in[0]);//out[0] = (2 * (in[0] + 2) * exp(-((in[0] + 2) * (in[0] + 2)) / 2)) + ((in[0] - 2) * exp(-((in[0] - 2) * (in[0] - 2)) / 2)) + ((in[0] + in[1]) / 10);
    out[1] = (2 * in[1]);//out[1] = (in[0] / 10) + ((5 * in[1] * in[1] * in[1] * in[1]) * exp(in[1] * in[1] * in[1] * in[1] * in[1]) - exp(-in[1] * in[1] * in[1] * in[1] * in[1]));
}

double fn2(double* in) {
 
    return (-2 * exp(-((in[0] + 2) * (in[0] + 2))/2)) - exp(-((in[0] - 2) * (in[0] - 2)) / 2) + ((in[0] * in[0]) / 20) + ((in[0] * in[1])/10) + exp(in[1] * in[1] * in[1] * in[1] * in [1]) + exp(-in[1] * in[1] * in[1] * in[1] * in[1]);
}

function<double(double*)> func = [](double* in) { return (-2 * exp(-((in[0] + 2) * (in[0] + 2)) / 2)) - exp(-((in[0] - 2) * (in[0] - 2)) / 2) + ((in[0] * in[0]) / 20) + ((in[0] * in[1]) / 10) + exp(in[1] * in[1] * in[1] * in[1] * in[1]) + exp(-in[1] * in[1] * in[1] * in[1] * in[1]);};
function<void(double*, double*)> grad = [](double* in, double* out) {
    out[0] = (2 * (in[0] + 2) * exp(-((in[0] + 2) * (in[0] + 2)) / 2)) + ((in[0] - 2) * exp(-((in[0] - 2) * (in[0] - 2)) / 2)) + ((in[0] + in[1]) / 10);
    out[1] = (in[0] / 10) + ((5 * in[1] * in[1] * in[1] * in[1]) * (exp(in[1] * in[1] * in[1] * in[1] * in[1]) - exp(-in[1] * in[1] * in[1] * in[1] * in[1])));
};
void fn_grad2(double* in, double* out) {

    out[0] = (2 * (in[0] + 2) * exp(-((in[0] + 2) * (in[0] + 2)) / 2)) + ((in[0] - 2) * exp(-((in[0] - 2) * (in[0] - 2)) / 2)) + ((in[0] + in[1]) / 10);
    out[1] = (in[0] / 10) + ((5 * in[1] * in[1] * in[1] * in[1]) * (exp(in[1] * in[1] * in[1] * in[1] * in[1]) - exp(-in[1] * in[1] * in[1] * in[1] * in[1])));
}



int main() {
    /*
    int const kNumOfVar = 2;
    double in[kNumOfVar];
    in[0] = 1;
    in[1] = 1;
    int m = 5;
    int DOF = kNumOfVar;
    double out[kNumOfVar];
    double step_size = 0.01;
    bool found_minimum = false;
    fn_grad(in, out);
    */
   
    /*
    SteepestDescent<double> steepestdescent(DOF);
    int i = 0;
    int maxcount = 1000;
    while (!found_minimum && i < maxcount) {
        ++i;
        steepestdescent.minimize_step(in, out, step_size);
        printf("f(%f, %f) = %f\n", in[0], in[1], fn(in));
        printf("grad of f(%f, %f) = [ %f, %f]\n", in[0], in[1], out[0], out[1]);
        fn_grad(in, out);
        double sum_of_residuals = 0;
        for (int i = 0; i < kNumOfVar; ++i) {
            sum_of_residuals += out[i] * out[i];

        }
        sum_of_residuals = sqrtf((sum_of_residuals / kNumOfVar));
        if (sum_of_residuals < 1e-5) {
            std::cout << "found minimum";
            found_minimum = true;
        }

    }*/
    
    /*
    LBFGS<double> lbfgs(m, DOF);
    lbfgs.init(in, out, step_size);
    fn_grad(in, out);
    int i = 0;
    int maxcount = 100000;
    while (!found_minimum && i < maxcount) {
        ++i;
        lbfgs.minimize_step(in, out, step_size);
        printf("f(%f, %f) = %f\n", in[0], in[1], fn(in));
        printf("grad of f(%f, %f) = [ %f, %f]\n", in[0], in[1], out[0], out[1]);
        fn_grad(in, out);
        double sum_of_residuals = 0;
        for (int i = 0; i < kNumOfVar; ++i) {
            sum_of_residuals += out[i] * out[i];

        }
        sum_of_residuals = sqrtf((sum_of_residuals / kNumOfVar));
        if (lbfgs.is_minimized() == true || sum_of_residuals < 1e-5) {
            std::cout << "found minimum\n\n";
            found_minimum = true;
        }

    }*/
   
    //second case

    int const kNumOfVar2 = 2;
    double in2[kNumOfVar2];
    in2[0] = 10;
    in2[1] = 1;
    int m2 = 5;
    int DOF2 = kNumOfVar2;
    double out2[kNumOfVar2];
    double step_size2 = 0.0001;
    bool found_minimum2 = false;
    fn_grad2(in2, out2);
    
    /*
    SteepestDescent<double> steepestdescent(DOF2);
    steepestdescent.init(in2, out2, step_size2);
    fn_grad2(in2, out2);
    printf("f(%f, %f) = %f\n", in2[0], in2[1], fn2(in2));
    printf("grad of f(%f, %f) = [ %f, %f]\n", in2[0], in2[1], out2[0], out2[1]);
    int i = 0;
    int maxcount2 = 100000;
    while (!found_minimum2 && i < maxcount2) {
        ++i;
        steepestdescent.minimize_step(in2, out2);
        printf("f(%f, %f) = %f\n", in2[0], in2[1], fn2(in2));
        fn_grad2(in2, out2);
        printf("grad of f(%f, %f) = [ %f, %f]\n", in2[0], in2[1], out2[0], out2[1]);
        double sum_of_residuals = 0;
        for (int i = 0; i < kNumOfVar2; ++i) {
            sum_of_residuals += out2[i] * out2[i];

        }
        sum_of_residuals = sqrtf((sum_of_residuals / kNumOfVar2));
        if (sum_of_residuals < 1e-5) {
            std::cout << "found minimum";
            found_minimum2 = true;
        }

    }
    */

   
    LBFGS<double> lbfgs2(m2, DOF2, func, grad);
    lbfgs2.init(in2, out2, step_size2);
    lbfgs2.minimize();
    /*
    fn_grad2(in2, out2);
    printf("f(%f, %f) = %f\n", in2[0], in2[1], fn2(in2));
    printf("grad of f(%f, %f) = [ %f, %f]\n", in2[0], in2[1], out2[0], out2[1]);
    int j = 0;
    int maxcount2 = 1000000;
    while (!found_minimum2 && j < maxcount2) {
        ++j;
        lbfgs2.minimize_step(in2, out2);
        printf("f(%f, %f) = %f\n", in2[0], in2[1], fn2(in2));
        fn_grad2(in2, out2);
        printf("grad of f(%f, %f) = [ %f, %f]\n", in2[0], in2[1], out2[0], out2[1]);
        double sum_of_residuals2 = 0;
        for (int i = 0; i < kNumOfVar2; ++i) {
            sum_of_residuals2 += out2[i] * out2[i];

        }
        sum_of_residuals2 = sqrtf((sum_of_residuals2 / kNumOfVar2));
        if (lbfgs2.is_minimized() == true || sum_of_residuals2 < 1e-7) {
            std::cout << "found minimum";
            found_minimum2 = true;
        }

    }
   */

 
    return 0;
   
}

