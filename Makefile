# Copyright 11/30/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine  #

#### Be sure that you have installed GMP with c++ interface ####

#### You may need to specify the include and lib paths of GMP ####
INC = # -I/usr/include
LIB = # -L/usr/lib

SHELL = /bin/sh
CC = g++
CFLAGS = -lm -lgmp -lgmpxx

TARGET = bach
SRC = $(shell echo src/*.cpp)
HEADER = $(shell echo src/*.h)
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(OBJ): %.o: %.cpp
	$(CC) -c -o $@ $< $(INC) 

$(TARGET): ${OBJ}
	$(CC) -o $@ $^ $(CFLAGS) $(INC) $(LIB)

.PHONY: clean

clean:
	-rm -f bach $(OBJ)

