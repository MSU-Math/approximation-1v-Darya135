#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>

class Window : public QWidget
{
    Q_OBJECT

  private:
    int func_id;
    const char *f_name;
    double a;
    double b;
    int n;
    int k;

    double (*f)(double);
    double (*df)(double);
    double (*d2f)(double);

    int display;
    int scale;
    int perturbation;
    double max_f;

    /* кэш построенных приближений */
    bool cache_valid;
    double cache_a, cache_b;
    int cache_n, cache_k, cache_perturbation, cache_scale;
    bool cache_11, cache_36;
    int cache_alloc_n;
    double *cache_X;
    double *cache_F;
    double *cache_A1;
    double *cache_A2;
    double *cache_extra;

    void free_cache();
    void build_cache(double a_scaled, double b_scaled);

    void DrawingFunction(QPainter &painter, double a, double b, double dx);
    void DrawingApproximation(QPainter &painter, double a, double b, double dx,
                              int n, const double *X, const double *F,
                              const double *A, int m);
    void DrawingError(QPainter &painter, double a, double b, double dx,
                      int n, const double *X, const double *F,
                      const double *A, int m);

  public:
    Window(QWidget *parent);
    ~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    int parse_command_line(int argc, char *argv[]);

  public slots:
    void change_func();

  protected:
    void paintEvent(QPaintEvent *event);
    void keyPressEvent(QKeyEvent *event);
};

#endif
