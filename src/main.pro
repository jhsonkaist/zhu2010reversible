
DEPENDPATH += . ../include/vcglib
INCLUDEPATH += . ../include/vcglib
CONFIG += console c++11
TEMPLATE = app

CONFIG -= app_bundle

CONFIG(debug, debug|release) {
    DESTDIR = ..
}

CONFIG(release, debug|release) {
    DESTDIR = ..
}

TARGET = main
SOURCES += main.cpp ../include/vcglib/wrap/ply/plylib.cpp