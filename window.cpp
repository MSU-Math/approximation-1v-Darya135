#include "window.h"
#include "inter_app.h"
#include <QKeyEvent>
#include <QPainter>
#include <cmath>
#include <cstdio>
#include <ctime>

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define DEFAULT_K 0
#define EPS 1e-16
#define MAXIMAL_N 10000000

static double GetMin(double y1, double y2, bool fp, int nes)
{
    if (nes == 1 && fp) {
        return y2;
    }
    return (y1 < y2) ? y1 : y2;
}

static double GetMax(double y1, double y2, bool fp, int nes)
{
    if (nes == 1 && fp) {
        return y2;
    }
    return (y1 > y2) ? y1 : y2;
}

static double f_0(double x)  { return 1.0 + 0.0 * x; }
static double f_1(double x)  { return x; }
static double f_2(double x)  { return pow(x, 2); }
static double f_3(double x)  { return pow(x, 3); }
static double f_4(double x)  { return pow(x, 4); }
static double f_5(double x)
{
    if (fabs(x) > pow(2, 9)) return x * 0.0;
    return exp(x);
}
static double f_6(double x) { return 1 / (25 * x * x + 1); }

static double df_0(double x)  { (void)x; return 0.0; }
static double df_1(double x)  { (void)x; return 1.0; }
static double df_2(double x)  { return 2.0 * x; }
static double df_3(double x)  { return 3.0 * x * x; }
static double df_4(double x)  { return 4.0 * x * x * x; }
static double df_5(double x)
{
    if (fabs(x) > pow(2, 9)) return 0.0;
    return exp(x);
}
static double df_6(double x)
{
    double d = 25.0 * x * x + 1.0;
    return -50.0 * x / (d * d);
}

static double d2f_0(double x) { (void)x; return 0.0; }
static double d2f_1(double x) { (void)x; return 0.0; }
static double d2f_2(double x) { (void)x; return 2.0; }
static double d2f_3(double x) { return 6.0 * x; }
static double d2f_4(double x) { return 12.0 * x * x; }
static double d2f_5(double x)
{
    if (fabs(x) > pow(2, 9)) return 0.0;
    return exp(x);
}
static double d2f_6(double x)
{
    double d = 25.0 * x * x + 1.0;
    return (3750.0 * x * x - 50.0) / (d * d * d);
}

Window::Window(QWidget *parent) : QWidget(parent)
{
    a = DEFAULT_A;
    b = DEFAULT_B;
    n = DEFAULT_N;
    k = DEFAULT_K;
    func_id = k;
    f = f_0;
    df = df_0;
    d2f = d2f_0;
    f_name = "f (x) = 1";
    display = 0;
    scale = 0;
    perturbation = 0;
    max_f = 0.0;

    cache_valid = false;
    cache_a = cache_b = 0.0;
    cache_n = cache_k = cache_perturbation = cache_scale = 0;
    cache_11 = cache_36 = false;
    cache_alloc_n = 0;
    cache_X = cache_F = cache_A1 = cache_A2 = cache_extra = 0;

    setFocusPolicy(Qt::StrongFocus);
}

Window::~Window()
{
    free_cache();
}

void Window::free_cache()
{
    delete[] cache_X;     cache_X = 0;
    delete[] cache_F;     cache_F = 0;
    delete[] cache_A1;    cache_A1 = 0;
    delete[] cache_A2;    cache_A2 = 0;
    delete[] cache_extra; cache_extra = 0;
    cache_alloc_n = 0;
    cache_valid = false;
}

void Window::build_cache(double a_scaled, double b_scaled)
{
    clock_t t_start = clock();

    /* Переиспользуем память, если n не изменилось.
     * Ключевое: не даём ядру заново размечать ~240 МБ
     * на каждом rebuild. */
    if (cache_alloc_n != n) {
        delete[] cache_X;     cache_X = new double[n];
        delete[] cache_F;     cache_F = new double[n];
        delete[] cache_A1;    cache_A1 = new double[n];
        delete[] cache_A2;    cache_A2 = new double[n];
        delete[] cache_extra; cache_extra = new double[n];
        cache_alloc_n = n;
    }

    double step = (b_scaled - a_scaled) / (n - 1);
    for (int i = 0; i < n; i++) {
        cache_X[i] = a_scaled + i * step;
        cache_F[i] = f(cache_X[i]);
        if (i == n / 2 && max_f > EPS) {
            cache_F[i] += perturbation * 0.1 * max_f;
        }
    }

    bool need_11 = (display == 0 || display == 2 || display == 3);
    bool need_36 = (display == 1 || display == 2 || display == 3);

    cache_11 = cache_36 = false;
    if (need_11 && n >= 2) {
        if (BuildingMethod11(n, cache_X, cache_F, cache_A1, cache_extra,
                             df(cache_X[0]), df(cache_X[n - 1])) == 0) {
            cache_11 = true;
        }
    }
    if (need_36 && n >= 2) {
        if (BuildingMethod36(n, cache_X, cache_F, cache_A2, cache_extra,
                             d2f(cache_X[0]), d2f(cache_X[n - 1])) == 0) {
            cache_36 = true;
        }
    }

    cache_a = a;
    cache_b = b;
    cache_n = n;
    cache_k = func_id;
    cache_perturbation = perturbation;
    cache_scale = scale;
    cache_valid = true;

    printf("build n=%d : %.3f s (CPU)\n", n,
           (double)(clock() - t_start) / CLOCKS_PER_SEC);
}

QSize Window::minimumSizeHint() const { return QSize(100, 100); }
QSize Window::sizeHint() const { return QSize(1000, 1000); }

int Window::parse_command_line(int argc, char *argv[])
{
    if (argc != 5) {
        qWarning("Wrong amount of arguments\n");
        return -1;
    }
    if (sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 ||
        b - a < 1.e-6) {
        qWarning("Wrong arguments a or b");
        return -2;
    }
    if (sscanf(argv[3], "%d", &n) != 1 || n < 2 || n > MAXIMAL_N) {
        qWarning("Wrong argument n");
        return -2;
    }
    if (sscanf(argv[4], "%d", &k) != 1 || k < 0 || k > 6) {
        qWarning("Wrong argument k");
        return -2;
    }
    func_id = k;
    change_func();
    return 0;
}

void Window::change_func()
{
    func_id = func_id % 7;

    switch (func_id) {
    case 0:
        f_name = "f (x) = 1";  f = f_0; df = df_0; d2f = d2f_0; break;
    case 1:
        f_name = "f (x) = x";  f = f_1; df = df_1; d2f = d2f_1; break;
    case 2:
        f_name = "f (x) = x^2"; f = f_2; df = df_2; d2f = d2f_2; break;
    case 3:
        f_name = "f (x) = x^3"; f = f_3; df = df_3; d2f = d2f_3; break;
    case 4:
        f_name = "f (x) = x^4"; f = f_4; df = df_4; d2f = d2f_4; break;
    case 5:
        f_name = "f (x) = e^x"; f = f_5; df = df_5; d2f = d2f_5; break;
    case 6:
        f_name = "f (x) = 1/(25*x^2 + 1)"; f = f_6; df = df_6; d2f = d2f_6; break;
    }
    cache_valid = false;
    update();
}

void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    double s = pow(2.0, scale);
    double a_scaled = a / s;
    double b_scaled = b / s;

    double delta_x = (b_scaled - a_scaled) / n;
    int draw_points = 1000;
    if (n <= 50) {
        draw_points = n * 10;
    }
    if (draw_points < 100) {
        draw_points = 100;
    }
    delta_x = (b_scaled - a_scaled) / draw_points;

    /* max|f| на видимом интервале */
    double max_f_global = 0.0;
    for (double xx = a_scaled; xx <= b_scaled; xx += delta_x) {
        double vv = fabs(f(xx));
        if (vv > max_f_global) {
            max_f_global = vv;
        }
    }
    max_f = max_f_global;
    printf("max|f| = %.16e\n", max_f);

    /* перестроить кэш только если параметры изменились */
    bool params_changed =
        !cache_valid || cache_a != a || cache_b != b || cache_n != n ||
        cache_k != func_id || cache_perturbation != perturbation ||
        cache_scale != scale;

    if (params_changed) {
        build_cache(a_scaled, b_scaled);
    }

    double *X = cache_X;
    double *F = cache_F;
    double *A1 = cache_A1;
    double *A2 = cache_A2;
    bool is_first_method = cache_11;
    bool is_second_method = cache_36;

    double max_abs = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
    double value = 0.0;
    double error = 0.0;
    double temp = 0.0;
    double error1 = 0.0;
    double error2 = 0.0;
    bool first_p = true;

    for (double x = a_scaled; x <= b_scaled; x += delta_x) {
        switch (display) {
        case 0:
            value = f(x);
            min_y = GetMin(min_y, value, first_p, 1);
            max_y = GetMax(max_y, value, first_p, 1);
            if (is_first_method) {
                value = EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1);
                min_y = GetMin(min_y, value, first_p, 0);
                max_y = GetMax(max_y, value, first_p, 0);
            }
            break;
        case 1:
            value = f(x);
            min_y = GetMin(min_y, value, first_p, 1);
            max_y = GetMax(max_y, value, first_p, 1);
            if (is_second_method) {
                value = EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2);
                min_y = GetMin(min_y, value, first_p, 0);
                max_y = GetMax(max_y, value, first_p, 0);
            }
            break;
        case 2:
            value = f(x);
            min_y = GetMin(min_y, value, first_p, 1);
            max_y = GetMax(max_y, value, first_p, 1);
            if (is_first_method) {
                value = EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1);
                min_y = GetMin(min_y, value, first_p, 0);
                max_y = GetMax(max_y, value, first_p, 0);
            }
            if (is_second_method) {
                value = EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2);
                min_y = GetMin(min_y, value, first_p, 0);
                max_y = GetMax(max_y, value, first_p, 0);
            }
            break;
        case 3:
            if (is_first_method) {
                value = fabs(EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1) - f(x));
                min_y = GetMin(min_y, value, first_p, 1);
                max_y = GetMax(max_y, value, first_p, 1);
                first_p = false;
            }
            if (is_second_method) {
                value = fabs(EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2) - f(x));
                min_y = GetMin(min_y, value, first_p, 1);
                max_y = GetMax(max_y, value, first_p, 1);
                first_p = false;
            }
            break;
        }
        first_p = false;
    }

    for (int kk = 0; kk < 2; kk++) {
        double x = (kk == 0) ? a_scaled : b_scaled;
        if (display != 3) {
            value = f(x);
            min_y = GetMin(min_y, value, first_p, 0);
            max_y = GetMax(max_y, value, first_p, 0);
        }
        if (is_first_method) {
            if (display != 3) {
                value = EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1);
            } else {
                value = fabs(EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1) - f(x));
            }
            min_y = GetMin(min_y, value, first_p, 0);
            max_y = GetMax(max_y, value, first_p, 0);
        }
        if (is_second_method) {
            if (display != 3) {
                value = EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2);
            } else {
                value = fabs(EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2) - f(x));
            }
            min_y = GetMin(min_y, value, first_p, 0);
            max_y = GetMax(max_y, value, first_p, 0);
        }
    }

    max_abs = fabs(min_y);
    if (fabs(max_y) > max_abs) {
        max_abs = fabs(max_y);
    }

    double delta_y = 0.05 * (max_y - min_y);
    if (delta_y < EPS) {
        delta_y = 1.0;
    }
    min_y -= delta_y;
    max_y += delta_y;

    painter.save();
    painter.translate(0, height());
    painter.scale(width() / (b_scaled - a_scaled), -height() / (max_y - min_y));
    painter.translate(-a_scaled, -min_y);
    QPen pen("black");
    pen.setWidth(0);
    painter.setPen(pen);

    switch (display) {
    case 0:
        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
        if (is_first_method) {
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, F, A1, 1);
        }
        break;
    case 1:
        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
        if (is_second_method) {
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, F, A2, 2);
        }
        break;
    case 2:
        DrawingFunction(painter, a_scaled, b_scaled, delta_x);
        if (is_first_method) {
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, F, A1, 1);
        }
        if (is_second_method) {
            DrawingApproximation(painter, a_scaled, b_scaled, delta_x, n, X, F, A2, 2);
        }
        break;
    case 3:
        if (is_first_method) {
            DrawingError(painter, a_scaled, b_scaled, delta_x, n, X, F, A1, 1);
        }
        if (is_second_method) {
            DrawingError(painter, a_scaled, b_scaled, delta_x, n, X, F, A2, 2);
        }
        break;
    }
    pen.setColor("black");
    painter.setPen(pen);
    painter.drawLine(a_scaled, 0, b_scaled, 0);
    painter.drawLine(0, min_y, 0, max_y);
    painter.restore();

    painter.setPen("black");
    painter.drawText(10, 20, QString("k=%1, %2").arg(k).arg(f_name));
    painter.drawText(10, 40, QString("max|f| = %1").arg(max_f));
    painter.drawText(10, 60, QString("scale = %1").arg(scale));
    painter.drawText(10, 80, QString("n = %1").arg(n));
    painter.drawText(10, 100, QString("perturbation p = %1").arg(perturbation));
    painter.drawText(10, 120, QString("display = %1").arg(display));

    int pos = 160;
    if (display == 3) {
        if (is_first_method) {
            for (double x = a_scaled; x <= b_scaled; x += delta_x) {
                temp = fabs(EvaluationMethod11(x, a_scaled, b_scaled, n, X, F, A1) - f(x));
                if (temp > error1) error1 = temp;
            }
            painter.drawText(10, pos, QString("error (method 11) = %1").arg(error1));
            pos += 20;
        }
        if (is_second_method) {
            for (double x = a_scaled; x <= b_scaled; x += delta_x) {
                temp = fabs(EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2) - f(x));
                if (temp > error2) error2 = temp;
            }
            painter.drawText(10, pos, QString("error (method 36) = %1").arg(error2));
            pos += 20;
        }
    } else {
        double residual = 0.0;
        if (is_second_method) {
            for (double x = a_scaled; x <= b_scaled; x += delta_x) {
                error = fabs(EvaluationMethod36(x, a_scaled, b_scaled, n, X, F, A2) - f(x));
                if (error > residual) residual = error;
            }
            painter.drawText(10, pos, QString("residual (method 36) = %1").arg(residual));
            pos += 20;
        }
    }

    painter.setPen(Qt::darkGreen);
    painter.drawText(10, pos, "green: f(x)"); pos += 20;
    if (is_first_method) {
        painter.setPen(Qt::red);
        painter.drawText(10, pos, "red: method 11 (clamped)"); pos += 20;
    }
    if (is_second_method) {
        painter.setPen(Qt::blue);
        painter.drawText(10, pos, "blue: method 36 (natural)"); pos += 20;
    }
    painter.setPen(Qt::black);
}

void Window::DrawingFunction(QPainter &painter, double a, double b, double dx)
{
    double x1 = a;
    double y1 = f(x1);
    QPen pen("green");
    pen.setWidth(0);
    painter.setPen(pen);
    for (double x2 = x1 + dx; x2 <= b + dx * 0.5; x2 += dx) {
        painter.drawLine(QPointF(x1, y1), QPointF(x2, f(x2)));
        x1 = x2;
        y1 = f(x2);
    }
}

void Window::DrawingApproximation(QPainter &painter, double a, double b,
                                  double dx, int n, const double *X,
                                  const double *F, const double *A, int m)
{
    double x1 = a;
    double y1 = 0.0;
    double y2 = 0.0;
    QPen pen("black");
    if (m == 1) {
        pen = QColor("red");
    } else {
        pen = QColor("blue");
    }
    pen.setWidth(0);
    painter.setPen(pen);
    if (m == 1) {
        y1 = EvaluationMethod11(x1, a, b, n, X, F, A);
    } else {
        y1 = EvaluationMethod36(x1, a, b, n, X, F, A);
    }
    for (double x2 = x1 + dx; x2 <= b; x2 += dx) {
        if (m == 2) {
            y2 = EvaluationMethod36(x2, a, b, n, X, F, A);
        } else {
            y2 = EvaluationMethod11(x2, a, b, n, X, F, A);
        }
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2;
        y1 = y2;
    }
}

void Window::DrawingError(QPainter &painter, double a, double b, double dx,
                          int n, const double *X, const double *F,
                          const double *A, int m)
{
    double x1 = a;
    double y2 = 0.0;
    double y1 = 0.0;
    QPen pen("black");
    if (m == 1) {
        pen = QColor("red");
    } else {
        pen = QColor("blue");
    }
    pen.setWidth(0);
    painter.setPen(pen);
    if (m == 2) {
        y1 = fabs(EvaluationMethod36(x1, a, b, n, X, F, A) - f(x1));
    } else {
        y1 = fabs(EvaluationMethod11(x1, a, b, n, X, F, A) - f(x1));
    }
    for (double x2 = x1 + dx; x2 <= b; x2 += dx) {
        if (m == 2) {
            y2 = fabs(EvaluationMethod36(x2, a, b, n, X, F, A) - f(x2));
        } else {
            y2 = fabs(EvaluationMethod11(x2, a, b, n, X, F, A) - f(x2));
        }
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        x1 = x2;
        y1 = y2;
    }
}

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        k = (k + 1) % 7;
        func_id = k;
        change_func();
        break;
    case Qt::Key_1:
        display = (display + 1) % 4;
        cache_valid = false;
        update();
        break;
    case Qt::Key_2:
        scale++;
        if (scale > 20) scale = 20;
        cache_valid = false;
        update();
        break;
    case Qt::Key_3:
        scale--;
        if (scale < -10) scale = -10;
        cache_valid = false;
        update();
        break;
    case Qt::Key_4:
        if (n <= MAXIMAL_N / 2) {
            n *= 2;
        } else {
            n = MAXIMAL_N;
        }
        cache_valid = false;
        update();
        break;
    case Qt::Key_5:
        n /= 2;
        if (n < 2) n = 2;
        cache_valid = false;
        update();
        break;
    case Qt::Key_6:
        perturbation++;
        if (perturbation > 50) perturbation = 50;
        cache_valid = false;
        update();
        break;
    case Qt::Key_7:
        perturbation--;
        if (perturbation < -50) perturbation = -50;
        cache_valid = false;
        update();
        break;
    default:
        QWidget::keyPressEvent(event);
    }
}
