#-------------------------------------------------
#
# Project created by QtCreator 2016-07-22T15:12:08
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BCQT
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    betweenness.cpp \
    graph.cpp

HEADERS  += mainwindow.h \
    betweenness.h \
    graph.h \
    define.h

FORMS    += mainwindow.ui
macx {
    QMAKE_CXXFLAGS += -std=c++11

    _BOOST_PATH = /usr/local/Cellar/boost/1.60.0_1
    INCLUDEPATH += "$${_BOOST_PATH}/include/"
    LIBS += -L$${_BOOST_PATH}/lib
    ## Use only one of these:
    LIBS += -lboost_chrono-mt -lboost_system # using dynamic lib (not sure if you need that "-mt" at the end or not)
    #LIBS += $${_BOOST_PATH}/lib/libboost_chrono-mt.a # using static lib
}
