#include "approximation.h"

#include <cmath>

static double node_coordinate(double left, double h, int index)
{
    return left + h * (double)index;
}

static double node_value(int index, double left, double h, int n, int func_id,
                         int perturbation, double perturbation_step,
                         FunctionPointer function_value)
{
    int mid;
    double x;
    double value;

    mid = n / 2;
    x = node_coordinate(left, h, index);
    value = function_value(func_id, x);

    if (index == mid) {
        value += (double)perturbation * perturbation_step;
    }

    return value;
}

static int interval_index(double point, double left, double right, int n)
{
    int index;
    double relative_index;

    if (point <= left) {
        return 0;
    }
    if (point >= right) {
        return n - 2;
    }

    relative_index = (point - left) * (double)(n - 1) / (right - left);
    index = (int)relative_index;

    if (index < 0) {
        index = 0;
    }
    if (index > n - 2) {
        index = n - 2;
    }

    return index;
}

static double evaluate_cubic(double point, double x0, double x1, double f0,
                             double f1, double d0, double d1)
{
    double h;
    double divided_difference;
    double c0;
    double c1;
    double c2;
    double c3;
    double dx;

    h = x1 - x0;
    divided_difference = (f1 - f0) / h;

    c0 = f0;
    c1 = d0;
    c2 = (3.0 * divided_difference - 2.0 * d0 - d1) / h;
    c3 = (d0 + d1 - 2.0 * divided_difference) / (h * h);

    dx = point - x0;
    return ((c3 * dx + c2) * dx + c1) * dx + c0;
}

int build_method_11(int n, double left, double right, int func_id,
                    int perturbation, double perturbation_step,
                    FunctionPointer function_value,
                    FunctionPointer function_first_derivative,
                    double *derivatives, double *temporary)
{
    int i;
    int row_count;
    double h;
    double eps;
    double value_im1;
    double value_i;
    double value_ip1;
    double delta_left;
    double delta_right;
    double rhs_value;
    double denominator;
    double lower;
    double upper;
    double x_left;
    double x_right;

    if (n < 2 || !(right > left) || function_value == 0 ||
        function_first_derivative == 0 || derivatives == 0 || temporary == 0) {
        return APPROXIMATION_ERROR_INVALID_ARGUMENTS;
    }

    h = (right - left) / (double)(n - 1);
    eps = 1.0e-300;
    x_left = left;
    x_right = right;

    derivatives[0] = function_first_derivative(func_id, x_left);
    derivatives[n - 1] = function_first_derivative(func_id, x_right);

    if (n == 2) {
        return APPROXIMATION_SUCCESS;
    }

    row_count = n - 2;

    for (i = 0; i < row_count; ++i) {
        value_im1 = node_value(i, left, h, n, func_id, perturbation,
                               perturbation_step, function_value);
        value_i = node_value(i + 1, left, h, n, func_id, perturbation,
                             perturbation_step, function_value);
        value_ip1 = node_value(i + 2, left, h, n, func_id, perturbation,
                               perturbation_step, function_value);

        delta_left = (value_i - value_im1) / h;
        delta_right = (value_ip1 - value_i) / h;
        rhs_value = 3.0 * (delta_left + delta_right);

        lower = 1.0;
        upper = 1.0;

        if (i == 0) {
            rhs_value -= derivatives[0];
            lower = 0.0;
        }
        if (i == row_count - 1) {
            rhs_value -= derivatives[n - 1];
            upper = 0.0;
        }

        if (i == 0) {
            denominator = 4.0;
        } else {
            denominator = 4.0 - lower * temporary[i - 1];
        }

        if (std::fabs(denominator) < eps) {
            return APPROXIMATION_ERROR_ZERO_PIVOT;
        }

        temporary[i] = upper / denominator;
        if (i == 0) {
            derivatives[i + 1] = rhs_value / denominator;
        } else {
            derivatives[i + 1] =
                (rhs_value - lower * derivatives[i]) / denominator;
        }
    }

    for (i = row_count - 2; i >= 0; --i) {
        derivatives[i + 1] -= temporary[i] * derivatives[i + 2];
    }

    return APPROXIMATION_SUCCESS;
}

int build_method_36(int n, double left, double right, int func_id,
                    int perturbation, double perturbation_step,
                    FunctionPointer function_value,
                    FunctionPointer function_second_derivative,
                    double *derivatives, double *temporary)
{
    int i;
    double h;
    double eps;
    double x_left;
    double x_right;
    double value_im1;
    double value_i;
    double value_ip1;
    double value_nm2;
    double value_nm1;
    double delta_left;
    double delta_right;
    double rhs_value;
    double denominator;
    double ddf_left;
    double ddf_right;

    if (n < 2 || !(right > left) || function_value == 0 ||
        function_second_derivative == 0 || derivatives == 0 || temporary == 0) {
        return APPROXIMATION_ERROR_INVALID_ARGUMENTS;
    }

    h = (right - left) / (double)(n - 1);
    eps = 1.0e-300;
    x_left = left;
    x_right = right;

    value_im1 = node_value(0, left, h, n, func_id, perturbation,
                           perturbation_step, function_value);
    value_i = node_value(1, left, h, n, func_id, perturbation,
                         perturbation_step, function_value);

    ddf_left = function_second_derivative(func_id, x_left);
    delta_right = (value_i - value_im1) / h;
    rhs_value = 3.0 * delta_right - 0.5 * ddf_left * h;

    denominator = 2.0;
    if (std::fabs(denominator) < eps) {
        return APPROXIMATION_ERROR_ZERO_PIVOT;
    }
    temporary[0] = 1.0 / denominator;
    derivatives[0] = rhs_value / denominator;

    for (i = 1; i < n - 1; ++i) {
        value_ip1 = node_value(i + 1, left, h, n, func_id, perturbation,
                               perturbation_step, function_value);
        delta_left = (value_i - value_im1) / h;
        delta_right = (value_ip1 - value_i) / h;
        rhs_value = 3.0 * (delta_left + delta_right);

        denominator = 4.0 - temporary[i - 1];
        if (std::fabs(denominator) < eps) {
            return APPROXIMATION_ERROR_ZERO_PIVOT;
        }
        temporary[i] = 1.0 / denominator;
        derivatives[i] = (rhs_value - derivatives[i - 1]) / denominator;

        value_im1 = value_i;
        value_i = value_ip1;
    }

    value_nm2 = value_im1;
    value_nm1 = value_i;
    ddf_right = function_second_derivative(func_id, x_right);
    delta_left = (value_nm1 - value_nm2) / h;
    rhs_value = 3.0 * delta_left + 0.5 * ddf_right * h;

    denominator = 2.0 - temporary[n - 2];
    if (std::fabs(denominator) < eps) {
        return APPROXIMATION_ERROR_ZERO_PIVOT;
    }
    temporary[n - 1] = 0.0;
    derivatives[n - 1] = (rhs_value - derivatives[n - 2]) / denominator;

    for (i = n - 2; i >= 0; --i) {
        derivatives[i] -= temporary[i] * derivatives[i + 1];
    }

    return APPROXIMATION_SUCCESS;
}

double evaluate_method_11(double point, double left, double right, int n,
                          int func_id, int perturbation,
                          double perturbation_step,
                          FunctionPointer function_value,
                          const double *derivatives)
{
    int index;
    double h;
    double x0;
    double x1;
    double f0;
    double f1;
    double d0;
    double d1;

    if (n < 2 || !(right > left) || function_value == 0 || derivatives == 0) {
        return NAN;
    }

    h = (right - left) / (double)(n - 1);
    index = interval_index(point, left, right, n);
    x0 = node_coordinate(left, h, index);
    x1 = x0 + h;

    f0 = node_value(index, left, h, n, func_id, perturbation, perturbation_step,
                    function_value);
    f1 = node_value(index + 1, left, h, n, func_id, perturbation,
                    perturbation_step, function_value);
    d0 = derivatives[index];
    d1 = derivatives[index + 1];

    return evaluate_cubic(point, x0, x1, f0, f1, d0, d1);
}

double evaluate_method_36(double point, double left, double right, int n,
                          int func_id, int perturbation,
                          double perturbation_step,
                          FunctionPointer function_value,
                          const double *derivatives)
{
    int index;
    double h;
    double x0;
    double x1;
    double f0;
    double f1;
    double d0;
    double d1;

    if (n < 2 || !(right > left) || function_value == 0 || derivatives == 0) {
        return NAN;
    }

    h = (right - left) / (double)(n - 1);
    index = interval_index(point, left, right, n);
    x0 = node_coordinate(left, h, index);
    x1 = x0 + h;

    f0 = node_value(index, left, h, n, func_id, perturbation, perturbation_step,
                    function_value);
    f1 = node_value(index + 1, left, h, n, func_id, perturbation,
                    perturbation_step, function_value);
    d0 = derivatives[index];
    d1 = derivatives[index + 1];

    return evaluate_cubic(point, x0, x1, f0, f1, d0, d1);
}
