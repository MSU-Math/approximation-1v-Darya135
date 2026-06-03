QMAKE_CXXFLAGS += -Werror -Wall -Wextra
HEADERS = window.h \
          functions.h \
          approximation.h
SOURCES = main.cpp \
          window.cpp \
          functions.cpp \
          approximation.cpp
QT += widgets
