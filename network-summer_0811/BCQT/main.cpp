#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Window w;
    w.setFixedSize(850,850);
    w.show();

    return a.exec();
}
