#include <QAction>
#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>

#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QMainWindow *window;
    QMenuBar *tool_bar;
    Window *graph_area;
    QAction *action;
    int error_code;

    window = new QMainWindow;
    tool_bar = new QMenuBar(window);
    graph_area = new Window(window);

    error_code = graph_area->parse_command_line(argc, argv);
    if (error_code != 0) {
        qWarning("Usage: ./basic_graph a b n k");
        qWarning("where a,b are double, n >= 2, k = 0,...,6");
        delete window;
        return -1;
    }

    action = tool_bar->addAction("&Function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));

    action =
        tool_bar->addAction("&Mode", graph_area, SLOT(change_graph_mode()));
    action->setShortcut(QString("1"));

    action = tool_bar->addAction("Zoom &in", graph_area, SLOT(zoom_in()));
    action->setShortcut(QString("2"));

    action = tool_bar->addAction("Zoom &out", graph_area, SLOT(zoom_out()));
    action->setShortcut(QString("3"));

    action = tool_bar->addAction("Increase &n", graph_area, SLOT(increase_n()));
    action->setShortcut(QString("4"));

    action = tool_bar->addAction("Decrease n", graph_area, SLOT(decrease_n()));
    action->setShortcut(QString("5"));

    action = tool_bar->addAction("Increase perturbation", graph_area,
                                 SLOT(increase_perturbation()));
    action->setShortcut(QString("6"));

    action = tool_bar->addAction("Decrease perturbation", graph_area,
                                 SLOT(decrease_perturbation()));
    action->setShortcut(QString("7"));

    action = tool_bar->addAction("E&xit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    tool_bar->setMaximumHeight(30);
    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Approximation: methods 11 and 36");
    window->show();
    graph_area->setFocus();

    return app.exec();
}
