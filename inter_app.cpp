#include "inter_app.h"
#include <cmath>
#include <stdio.h>

#define EPS_PIVOT 1.0e-300

int BuildingMethod11(int n, const double *x, const double *f, double *a,
                     double *extra, double df_left, double df_right)
{
    if (n < 2 || !x || !f || !a || !extra) {
        return -1;
    }

    double h = (x[n - 1] - x[0]) / (double)(n - 1);
    double *d = a;
    double *temp = extra;

    d[0] = df_left;
    d[n - 1] = df_right;
    if (n == 2) {
        return 0;
    }

    int row_count = n - 2;
    for (int i = 0; i < row_count; ++i) {
        int idx = i + 1;
        double delta_left = (f[idx] - f[idx - 1]) / h;
        double delta_right = (f[idx + 1] - f[idx]) / h;
        double rhs = 3.0 * (delta_left + delta_right);

        double lower = 1.0;
        double upper = 1.0;
        if (i == 0) {
            rhs -= d[0];
            lower = 0.0;
        }
        if (i == row_count - 1) {
            rhs -= d[n - 1];
            upper = 0.0;
        }

        double denom = (i == 0) ? 4.0 : (4.0 - lower * temp[i - 1]);
        if (fabs(denom) < EPS_PIVOT) {
            return -1;
        }
        temp[i] = upper / denom;
        if (i == 0) {
            d[idx] = rhs / denom;
        } else {
            d[idx] = (rhs - lower * d[idx - 1]) / denom;
        }
    }

    for (int i = row_count - 2; i >= 0; --i) {
        d[i + 1] -= temp[i] * d[i + 2];
    }

    return 0;
}

int BuildingMethod36(int n, const double *x, const double *f, double *a,
                     double *extra, double d2f_left, double d2f_right)
{
    if (n < 2 || !x || !f || !a || !extra) {
        return -1;
    }

    double h = (x[n - 1] - x[0]) / (double)(n - 1);
    double *d = a;
    double *temp = extra;

    double rhs = 3.0 * (f[1] - f[0]) / h - 0.5 * d2f_left * h;
    double denom = 2.0;
    if (fabs(denom) < EPS_PIVOT) {
        return -1;
    }
    temp[0] = 1.0 / denom;
    d[0] = rhs / denom;

    for (int i = 1; i < n - 1; ++i) {
        rhs = 3.0 * (f[i + 1] - f[i - 1]) / h;
        denom = 4.0 - temp[i - 1];
        if (fabs(denom) < EPS_PIVOT) {
            return -1;
        }
        temp[i] = 1.0 / denom;
        d[i] = (rhs - d[i - 1]) / denom;
    }

    rhs = 3.0 * (f[n - 1] - f[n - 2]) / h + 0.5 * d2f_right * h;
    denom = 2.0 - temp[n - 2];
    if (fabs(denom) < EPS_PIVOT) {
        return -1;
    }
    d[n - 1] = (rhs - d[n - 2]) / denom;

    for (int i = n - 2; i >= 0; --i) {
        d[i] -= temp[i] * d[i + 1];
    }

    return 0;
}

static int FindSegment(double x, int n, const double *X)
{
    if (x <= X[0]) {
        return 0;
    }
    if (x >= X[n - 1]) {
        return n - 2;
    }
    int left = 0;
    int right = n - 2;
    while (left < right) {
        int middle = (left + right) / 2;
        if (X[middle + 1] < x) {
            left = middle + 1;
        } else {
            right = middle;
        }
    }
    return left;
}

static double EvalCubic(int i, double x, const double *X, const double *F,
                        const double *d)
{
    double h = X[i + 1] - X[i];
    double dx = x - X[i];
    double delta = (F[i + 1] - F[i]) / h;

    double a0 = F[i];
    double a1 = d[i];
    double a2 = (3.0 * delta - 2.0 * d[i] - d[i + 1]) / h;
    double a3 = (d[i] + d[i + 1] - 2.0 * delta) / (h * h);

    return a0 + a1 * dx + a2 * dx * dx + a3 * dx * dx * dx;
}

double EvaluationMethod11(double x, double a, double b, int n,
                          const double *X, const double *F, const double *A)
{
    (void)a;
    (void)b;
    if (n < 2 || !X || !F || !A) {
        return 0.0;
    }
    int i = FindSegment(x, n, X);
    return EvalCubic(i, x, X, F, A);
}

double EvaluationMethod36(double x, double a, double b, int n,
                          const double *X, const double *F, const double *A)
{
    (void)a;
    (void)b;
    if (n < 2 || !X || !F || !A) {
        return 0.0;
    }
    int i = FindSegment(x, n, X);
    return EvalCubic(i, x, X, F, A);
}
