
DEPENDPATH += . ../include/vcglib
INCLUDEPATH += . ../include/vcglib
 
TEMPLATE = app

QT += core
QT -= gui

CONFIG += c++11
CONFIG += console
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++11

CONFIG(debug, debug|release) {
    DESTDIR = ..
}

CONFIG(release, debug|release) {
    DESTDIR = ..
}

TARGET = temp
SOURCES += main.cpp ../include/vcglib/wrap/ply/plylib.cpp