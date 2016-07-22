#include "mainwindow.h"

Window::Window()
  : QOpenGLWidget()
{
  QSurfaceFormat format;
  format.setSamples(8);
  format.setStencilBufferSize(8);
  setFormat(format);
  printf("Version Major:%d minor:%d \n",format.version().first, format.version().second);
  datachanged = colorchanged = true;
}

void Window::updateWindow(){
  datachanged = colorchanged = true;
  update();
}

void Window::updateData(){
  datachanged = true;
  update();
}

void Window::updateColor(){
  colorchanged = true;
  update();
}

void Window::initializeGL(){
    std::cout << "here"<<std::endl;
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClearDepthf(1.0f);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
}

void Window::resizeGL(int width, int height){
  (void) width;
  (void) height;
}

void Window::paintGL(){

}
