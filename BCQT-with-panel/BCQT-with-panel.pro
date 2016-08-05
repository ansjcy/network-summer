#-------------------------------------------------
#
# Project created by QtCreator 2016-07-25T14:44:00
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = BCQT-with-panel
TEMPLATE = app


SOURCES += main.cpp\
        widget.cpp \
    betweenness.cpp \
    graph.cpp

HEADERS  += widget.h \
    betweenness.h \
    graph.h \
    define.h

FORMS    += widget.ui

macx {

# remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE += -O3

    QMAKE_CXXFLAGS += -std=c++11

    _BOOST_PATH = /usr/local/Cellar/boost/1.60.0_1
    INCLUDEPATH += "$${_BOOST_PATH}/include/"
    LIBS += -L$${_BOOST_PATH}/lib
    ## Use only one of these:
    LIBS += -lboost_chrono-mt -lboost_system # using dynamic lib (not sure if you need that "-mt" at the end or not)
    #LIBS += $${_BOOST_PATH}/lib/libboost_chrono-mt.a # using static lib
}
