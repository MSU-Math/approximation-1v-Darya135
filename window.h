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

    double (*f)(double); /* сама функция */
    double (*df)(double); /* первая производная  — для метода 11 */
    double (*d2f)(double); /* вторая производная  — для метода 36 */

    int display; /* 0..3 — что показывать */
    int scale; /* масштаб по X */
    int perturbation; /* p — возмущение f(x_{n/2}) */
    double max_f; /* max|f| на [a, b] — для возмущения */

    void DrawingFunction(QPainter &painter, double a, double b, double dx);
    void DrawingApproximation(QPainter &painter, double a, double b, double dx,
                              int n, const double *X, const double *A, int m);
    void DrawingError(QPainter &painter, double a, double b, double dx, int n,
                      const double *X, const double *A, int m);

  public:
    explicit Window(QWidget *parent);

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
