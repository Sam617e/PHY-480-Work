TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -LC:/Armadillo/include/newblas/old -lblas_win64_MT
LIBS += -LC:/Armadillo/include/newblas/old -llapack_win64_MT

INCLUDEPATH += C:/Armadillo/include
DEPENDPATH += C:/Armadillo/include
SOURCES += \
    main.cpp
