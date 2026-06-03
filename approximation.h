#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#define APPROXIMATION_SUCCESS 0
#define APPROXIMATION_ERROR_INVALID_ARGUMENTS -1
#define APPROXIMATION_ERROR_GRID -2
#define APPROXIMATION_ERROR_ZERO_PIVOT -3

typedef double (*FunctionPointer)(int, double);

int build_method_11(int n, double left, double right, int func_id,
                    int perturbation, double perturbation_step,
                    FunctionPointer function_value,
                    FunctionPointer function_first_derivative,
                    double *derivatives, double *temporary);

int build_method_36(int n, double left, double right, int func_id,
                    int perturbation, double perturbation_step,
                    FunctionPointer function_value,
                    FunctionPointer function_second_derivative,
                    double *derivatives, double *temporary);

double evaluate_method_11(double point, double left, double right, int n,
                          int func_id, int perturbation,
                          double perturbation_step,
                          FunctionPointer function_value,
                          const double *derivatives);

double evaluate_method_36(double point, double left, double right, int n,
                          int func_id, int perturbation,
                          double perturbation_step,
                          FunctionPointer function_value,
                          const double *derivatives);

#endif
