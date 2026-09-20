#ifndef INTER_APP_H
#define INTER_APP_H

int BuildingMethod11(int n, const double *x, const double *f, double *a,
                     double *extra, double df_left, double df_right);

double EvaluationMethod11(double x, double a, double b, int n, const double *X,
                          const double *F, const double *A);

int BuildingMethod36(int n, const double *x, const double *f, double *a,
                     double *extra, double d2f_left, double d2f_right);

double EvaluationMethod36(double x, double a, double b, int n, const double *X,
                          const double *F, const double *A);

#endif
