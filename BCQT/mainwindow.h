#ifndef WINDOW_H
#define WINDOW_H
#include <QOpenGLWidget>
#include <glu.h>
#include <OpenGL.h>
#include <iostream>

class Window : public QOpenGLWidget {
  Q_OBJECT
public:
    Window();

    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int width, int height);

    void updateWindow();
    void updateData();
    void updateColor();

protected:
    bool datachanged, colorchanged;

   //add by me


signals:
    void emitUpdate();
    void emitColorUpdate();
    void emitDataUpdate();

};

#endif // WINDOW_H
