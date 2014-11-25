OBJS = main.o serial.o parallel.o
CC = g++
DEBUG = 
CFLAGS = $(DEBUG) -std=c++11
LIBS = -lgomp -fopenmp -lm -pthread

all: pearson

pearson: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o pearson $(LIBS)

main.o: main.cpp
	$(CC) -c $(CFLAGS) main.cpp $(LIBS)

serial.o: serial.cpp
	$(CC) -c $(CFLAGS) serial.cpp $(LIBS)

parallel.o: parallel.cpp
	$(CC) -c $(CFLAGS) parallel.cpp $(LIBS)

clean:
	rm -rf *o pearson