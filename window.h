#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

#define GRAPH_MODE_METHOD_11 0
#define GRAPH_MODE_METHOD_36 1
#define GRAPH_MODE_BOTH 2
#define GRAPH_MODE_ERRORS 3

class Window : public QWidget
{
    Q_OBJECT

  private:
    double a;
    double b;
    int n;
    int func_id;
    int graph_mode;
    int scale_power;
    int perturbation;
    int approximation_error;

    double *derivatives_11;
    double *derivatives_36;

    void free_memory();
    int allocate_memory();
    int rebuild_approximation();
    double exact_value(double x) const;
    double method_11_value(double x) const;
    double method_36_value(double x) const;
    double perturbation_step() const;
    void get_visible_interval(double *left, double *right) const;
    const char *graph_mode_name() const;

  public:
    explicit Window(QWidget *parent);
    ~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    int parse_command_line(int argc, char *argv[]);

  public slots:
    void change_func();
    void change_graph_mode();
    void zoom_in();
    void zoom_out();
    void increase_n();
    void decrease_n();
    void increase_perturbation();
    void decrease_perturbation();

  protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
};

#endif
