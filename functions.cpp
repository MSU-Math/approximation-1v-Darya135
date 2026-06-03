#include "functions.h"

#include <cmath>

static double max_double(double x, double y)
{
    if (x > y) {
        return x;
    }
    return y;
}

double function_value(int k, double x)
{
    double d;

    switch (k) {
    case 0:
        return 1.0;
    case 1:
        return x;
    case 2:
        return x * x;
    case 3:
        return x * x * x;
    case 4:
        return x * x * x * x;
    case 5:
        return std::exp(x);
    case 6:
        d = 25.0 * x * x + 1.0;
        return 1.0 / d;
    default:
        return 0.0;
    }
}

double function_first_derivative(int k, double x)
{
    double d;

    switch (k) {
    case 0:
        return 0.0;
    case 1:
        return 1.0;
    case 2:
        return 2.0 * x;
    case 3:
        return 3.0 * x * x;
    case 4:
        return 4.0 * x * x * x;
    case 5:
        return std::exp(x);
    case 6:
        d = 25.0 * x * x + 1.0;
        return -50.0 * x / (d * d);
    default:
        return 0.0;
    }
}

double function_second_derivative(int k, double x)
{
    double d;

    switch (k) {
    case 0:
        return 0.0;
    case 1:
        return 0.0;
    case 2:
        return 2.0;
    case 3:
        return 6.0 * x;
    case 4:
        return 12.0 * x * x;
    case 5:
        return std::exp(x);
    case 6:
        d = 25.0 * x * x + 1.0;
        return (3750.0 * x * x - 50.0) / (d * d * d);
    default:
        return 0.0;
    }
}

double function_max_abs(int k, double a, double b)
{
    double v1;
    double v2;

    v1 = std::fabs(function_value(k, a));
    v2 = std::fabs(function_value(k, b));

    if (k == 0) {
        return 1.0;
    }
    if (k == 2 || k == 4) {
        return max_double(v1, v2);
    }
    if (k == 5) {
        return max_double(v1, v2);
    }
    if (k == 6 && a <= 0.0 && b >= 0.0) {
        return 1.0;
    }

    return max_double(v1, v2);
}

const char *function_name(int k)
{
    switch (k) {
    case 0:
        return "1";
    case 1:
        return "x";
    case 2:
        return "x^2";
    case 3:
        return "x^3";
    case 4:
        return "x^4";
    case 5:
        return "exp(x)";
    case 6:
        return "1/(25*x^2+1)";
    default:
        return "unknown";
    }
}
