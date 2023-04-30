CC = g++
CFLAGS = -g -Wall -std=c++11
TARGET = blcbt

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp