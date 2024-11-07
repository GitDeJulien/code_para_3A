CPP = g++
CFLAGS = -Wall -Wextra -O3 -std=c++14

SRCS = source/data.cpp source/functions.cpp source/main.cpp
TARGET = run

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CPP) $(CFLAGS) -o source/$(TARGET) $(SRCS)


exec:
	./source/run "./input/data.dat"

clean:
	rm -f source/*.o
	rm -f output/*.dat

clear:
	rm -f $(TARGET)