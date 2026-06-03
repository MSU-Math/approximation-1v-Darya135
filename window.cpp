#include "window.h"

#include "approximation.h"
#include "functions.h"

#include <QKeyEvent>
#include <QPainter>
#include <QPolygonF>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>

#define MINIMAL_N 2
#define FUNCTION_COUNT 7
#define SAMPLE_COUNT_MIN 100
#define SAMPLE_COUNT_MAX 1200
#define SCALE_POWER_MIN -8
#define SCALE_POWER_MAX 8
#define ERROR_DISPLAY_EPS 1.0e-12
#define MAXIMAL_N 10000000

Window::Window(QWidget *parent) : QWidget(parent)
{
    a = 0.0;
    b = 0.0;
    n = 0;
    func_id = 0;
    graph_mode = GRAPH_MODE_BOTH;
    scale_power = 0;
    perturbation = 0;
    approximation_error = APPROXIMATION_ERROR_INVALID_ARGUMENTS;

    derivatives_11 = 0;
    derivatives_36 = 0;

    setFocusPolicy(Qt::StrongFocus);
}

Window::~Window()
{
    free_memory();
}

QSize Window::minimumSizeHint() const
{
    return QSize(100, 100);
}

QSize Window::sizeHint() const
{
    return QSize(1000, 800);
}

void Window::free_memory()
{
    delete[] derivatives_11;
    delete[] derivatives_36;
    derivatives_11 = 0;
    derivatives_36 = 0;
}

int Window::allocate_memory()
{
    free_memory();

    if (n < MINIMAL_N) {
        return APPROXIMATION_ERROR_INVALID_ARGUMENTS;
    }

    derivatives_11 = new (std::nothrow) double[n];
    derivatives_36 = new (std::nothrow) double[n];
    if (derivatives_11 == 0 || derivatives_36 == 0) {
        free_memory();
        return APPROXIMATION_ERROR_INVALID_ARGUMENTS;
    }

    return APPROXIMATION_SUCCESS;
}

double Window::perturbation_step() const
{
    return 0.1 * function_max_abs(func_id, a, b);
}

int Window::rebuild_approximation()
{
    int error_code;
    double *temporary;
    double step;

    error_code = allocate_memory();
    if (error_code != APPROXIMATION_SUCCESS) {
        approximation_error = error_code;
        update();
        return error_code;
    }

    temporary = new (std::nothrow) double[n];
    if (temporary == 0) {
        free_memory();
        approximation_error = APPROXIMATION_ERROR_INVALID_ARGUMENTS;
        update();
        return approximation_error;
    }

    step = perturbation_step();

    error_code =
        build_method_11(n, a, b, func_id, perturbation, step, function_value,
                        function_first_derivative, derivatives_11, temporary);
    if (error_code == APPROXIMATION_SUCCESS) {
        error_code = build_method_36(n, a, b, func_id, perturbation, step,
                                     function_value, function_second_derivative,
                                     derivatives_36, temporary);
    }

    delete[] temporary;

    if (error_code != APPROXIMATION_SUCCESS) {
        approximation_error = error_code;
        update();
        return error_code;
    }

    approximation_error = APPROXIMATION_SUCCESS;
    update();
    return APPROXIMATION_SUCCESS;
}

int Window::parse_command_line(int argc, char *argv[])
{
    int parsed;

    if (argc != 5) {
        return -1;
    }

    parsed = sscanf(argv[1], "%lf", &a);
    if (parsed != 1) {
        return -1;
    }

    parsed = sscanf(argv[2], "%lf", &b);
    if (parsed != 1) {
        return -1;
    }

    parsed = sscanf(argv[3], "%d", &n);
    if (parsed != 1) {
        return -1;
    }

    parsed = sscanf(argv[4], "%d", &func_id);
    if (parsed != 1) {
        return -1;
    }

    if (!(b > a) || n < MINIMAL_N || n > MAXIMAL_N || func_id < 0 ||
        func_id >= FUNCTION_COUNT) {
        return -1;
    }

    return rebuild_approximation();
}

void Window::change_func()
{
    func_id = (func_id + 1) % FUNCTION_COUNT;
    perturbation = 0;
    scale_power = 0;
    rebuild_approximation();
}

void Window::change_graph_mode()
{
    graph_mode = (graph_mode + 1) % 4;
    update();
}

void Window::zoom_in()
{
    if (scale_power < SCALE_POWER_MAX) {
        ++scale_power;
        update();
    }
}

void Window::zoom_out()
{
    if (scale_power > SCALE_POWER_MIN) {
        --scale_power;
        update();
    }
}

void Window::increase_n()
{
    if (n <= MAXIMAL_N / 2) {
        n *= 2;
    } else {
        n = MAXIMAL_N;
    }
    rebuild_approximation();
}

void Window::decrease_n()
{
    if (n > MINIMAL_N) {
        n = (n + 1) / 2;
        if (n < MINIMAL_N) {
            n = MINIMAL_N;
        }
        rebuild_approximation();
    }
}

void Window::increase_perturbation()
{
    ++perturbation;
    rebuild_approximation();
}

void Window::decrease_perturbation()
{
    --perturbation;
    rebuild_approximation();
}

void Window::get_visible_interval(double *left, double *right) const
{
    double multiplier;
    double center;
    double half_length;

    multiplier = std::pow(2.0, (double)scale_power);
    center = 0.5 * (a + b);
    half_length = 0.5 * (b - a) / multiplier;

    *left = center - half_length;
    *right = center + half_length;
}

double Window::exact_value(double x) const
{
    return function_value(func_id, x);
}

double Window::method_11_value(double x) const
{
    return evaluate_method_11(x, a, b, n, func_id, perturbation,
                              perturbation_step(), function_value,
                              derivatives_11);
}

double Window::method_36_value(double x) const
{
    return evaluate_method_36(x, a, b, n, func_id, perturbation,
                              perturbation_step(), function_value,
                              derivatives_36);
}

const char *Window::graph_mode_name() const
{
    switch (graph_mode) {
    case GRAPH_MODE_METHOD_11:
        return "f and method 11";
    case GRAPH_MODE_METHOD_36:
        return "f and method 36";
    case GRAPH_MODE_BOTH:
        return "f and methods 11, 36";
    case GRAPH_MODE_ERRORS:
        return "errors of methods 11, 36";
    default:
        return "unknown";
    }
}

static double display_error_value(double value)
{
    if (std::fabs(value) < ERROR_DISPLAY_EPS) {
        return 0.0;
    }
    return value;
}

static bool is_inside_interval(double x, double left, double right)
{
    return x >= left && x <= right;
}

void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    QPen pen;
    QString info;
    int samples;
    int i;
    int text_y;
    double left;
    double right;
    double x;
    double y;
    double y_function;
    double y_method_11;
    double y_method_36;
    double min_y;
    double max_y;
    double max_abs;
    double y_current;
    bool initialized;
    QPolygonF points;

    painter.fillRect(rect(), Qt::white);

    if (approximation_error != APPROXIMATION_SUCCESS) {
        painter.setPen(Qt::red);
        painter.drawText(
            10, 25,
            QString("approximation error: %1").arg(approximation_error));
        return;
    }

    get_visible_interval(&left, &right);
    samples = width();
    if (samples < SAMPLE_COUNT_MIN) {
        samples = SAMPLE_COUNT_MIN;
    }
    if (samples > SAMPLE_COUNT_MAX) {
        samples = SAMPLE_COUNT_MAX;
    }

    initialized = false;
    min_y = 0.0;
    max_y = 0.0;

    for (i = 0; i <= samples; ++i) {
        x = left + (right - left) * (double)i / (double)samples;
        y_function = 0.0;
        y_method_11 = 0.0;
        y_method_36 = 0.0;

        if (graph_mode != GRAPH_MODE_ERRORS) {
            y_function = exact_value(x);
            if (!initialized) {
                min_y = y_function;
                max_y = y_function;
                initialized = true;
            }
            if (y_function < min_y) {
                min_y = y_function;
            }
            if (y_function > max_y) {
                max_y = y_function;
            }
        }

        if ((graph_mode == GRAPH_MODE_METHOD_11 ||
             graph_mode == GRAPH_MODE_BOTH) &&
            is_inside_interval(x, a, b)) {
            y_method_11 = method_11_value(x);
            if (y_method_11 < min_y) {
                min_y = y_method_11;
            }
            if (y_method_11 > max_y) {
                max_y = y_method_11;
            }
        }

        if ((graph_mode == GRAPH_MODE_METHOD_36 ||
             graph_mode == GRAPH_MODE_BOTH) &&
            is_inside_interval(x, a, b)) {
            y_method_36 = method_36_value(x);
            if (y_method_36 < min_y) {
                min_y = y_method_36;
            }
            if (y_method_36 > max_y) {
                max_y = y_method_36;
            }
        }

        if (graph_mode == GRAPH_MODE_ERRORS && is_inside_interval(x, a, b)) {
            y_function = exact_value(x);
            y_method_11 = method_11_value(x);
            y = display_error_value(y_method_11 - y_function);
            if (!initialized) {
                min_y = y;
                max_y = y;
                initialized = true;
            }
            if (y < min_y) {
                min_y = y;
            }
            if (y > max_y) {
                max_y = y;
            }

            y_method_36 = method_36_value(x);
            y = display_error_value(y_method_36 - y_function);
            if (y < min_y) {
                min_y = y;
            }
            if (y > max_y) {
                max_y = y;
            }
        }
    }

    max_abs = std::fabs(min_y);
    if (std::fabs(max_y) > max_abs) {
        max_abs = std::fabs(max_y);
    }

    if (graph_mode == GRAPH_MODE_ERRORS && max_abs < ERROR_DISPLAY_EPS) {
        min_y = -ERROR_DISPLAY_EPS;
        max_y = ERROR_DISPLAY_EPS;
    }

    if (min_y > 0.0) {
        min_y = 0.0;
    }
    if (max_y < 0.0) {
        max_y = 0.0;
    }

    if (min_y == max_y) {
        min_y -= 1.0;
        max_y += 1.0;
    }

    y = 0.05 * (max_y - min_y);
    min_y -= y;
    max_y += y;

    std::printf("mode=%s k=%d n=%d scale=%d perturbation=%d max_abs=%.16e\n",
                graph_mode_name(), func_id, n, scale_power, perturbation,
                max_abs);

    painter.save();
    painter.translate(0.5 * width(), 0.5 * height());
    painter.scale(width() / (right - left), -height() / (max_y - min_y));
    painter.translate(-0.5 * (left + right), -0.5 * (min_y + max_y));

    pen.setWidth(0);

    pen.setColor(Qt::red);
    painter.setPen(pen);
    painter.drawLine(QPointF(left, 0.0), QPointF(right, 0.0));
    painter.drawLine(QPointF(0.0, min_y), QPointF(0.0, max_y));

    if (graph_mode != GRAPH_MODE_ERRORS) {
        points.clear();
        points.reserve(samples + 1);
        pen.setColor(Qt::black);
        painter.setPen(pen);
        for (i = 0; i <= samples; ++i) {
            x = left + (right - left) * (double)i / (double)samples;
            y_current = exact_value(x);
            points.append(QPointF(x, y_current));
        }
        painter.drawPolyline(points);
    }

    if (graph_mode == GRAPH_MODE_METHOD_11 || graph_mode == GRAPH_MODE_BOTH ||
        graph_mode == GRAPH_MODE_ERRORS) {
        points.clear();
        points.reserve(samples + 1);
        pen.setColor(Qt::blue);
        painter.setPen(pen);
        for (i = 0; i <= samples; ++i) {
            x = left + (right - left) * (double)i / (double)samples;
            if (is_inside_interval(x, a, b)) {
                if (graph_mode == GRAPH_MODE_ERRORS) {
                    y_current = display_error_value(method_11_value(x) -
                                                    exact_value(x));
                } else {
                    y_current = method_11_value(x);
                }
                points.append(QPointF(x, y_current));
            } else if (points.size() > 1) {
                painter.drawPolyline(points);
                points.clear();
            } else {
                points.clear();
            }
        }
        if (points.size() > 1) {
            painter.drawPolyline(points);
        }
    }

    if (graph_mode == GRAPH_MODE_METHOD_36 || graph_mode == GRAPH_MODE_BOTH ||
        graph_mode == GRAPH_MODE_ERRORS) {
        points.clear();
        points.reserve(samples + 1);
        pen.setColor(Qt::darkGreen);
        painter.setPen(pen);
        for (i = 0; i <= samples; ++i) {
            x = left + (right - left) * (double)i / (double)samples;
            if (is_inside_interval(x, a, b)) {
                if (graph_mode == GRAPH_MODE_ERRORS) {
                    y_current = display_error_value(method_36_value(x) -
                                                    exact_value(x));
                } else {
                    y_current = method_36_value(x);
                }
                points.append(QPointF(x, y_current));
            } else if (points.size() > 1) {
                painter.drawPolyline(points);
                points.clear();
            } else {
                points.clear();
            }
        }
        if (points.size() > 1) {
            painter.drawPolyline(points);
        }
    }

    painter.restore();

    painter.setPen(Qt::black);
    text_y = 20;
    info = QString("k=%1  f(x)=%2").arg(func_id).arg(function_name(func_id));
    painter.drawText(10, text_y, info);
    text_y += 20;
    info = QString("n=%1  mode=%2").arg(n).arg(graph_mode_name());
    painter.drawText(10, text_y, info);
    text_y += 20;
    info = QString("scale=%1  interval=[%2,%3]")
               .arg(scale_power)
               .arg(left, 0, 'g', 4)
               .arg(right, 0, 'g', 4);
    painter.drawText(10, text_y, info);
    text_y += 20;
    info = QString("perturbation=%1  max_abs=%2")
               .arg(perturbation)
               .arg(max_abs, 0, 'e', 3);
    painter.drawText(10, text_y, info);
    text_y += 20;
    painter.drawText(10, text_y,
                     "keys: 0 func, 1 mode, 2/3 zoom, 4/5 n, 6/7 perturb");
    text_y += 20;
    painter.setPen(Qt::black);
    painter.drawText(10, text_y, "black: f(x)");
    text_y += 20;
    painter.setPen(Qt::blue);
    painter.drawText(10, text_y, "blue: method 11");
    text_y += 20;
    painter.setPen(Qt::darkGreen);
    painter.drawText(10, text_y, "green: method 36");
}

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        change_func();
        break;
    case Qt::Key_1:
        change_graph_mode();
        break;
    case Qt::Key_2:
        zoom_in();
        break;
    case Qt::Key_3:
        zoom_out();
        break;
    case Qt::Key_4:
        increase_n();
        break;
    case Qt::Key_5:
        decrease_n();
        break;
    case Qt::Key_6:
        increase_perturbation();
        break;
    case Qt::Key_7:
        decrease_perturbation();
        break;
    default:
        QWidget::keyPressEvent(event);
        break;
    }
}
