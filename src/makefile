EXEC = visualization
OBJS = main.o cell.o weather_data.o lodepng.o
SRCS = .
INCS = ./include
LIBS = ./include/lodepng

CC = g++
CFLAGS = -Wall -O3 -c -std=c++0x
LFLAGS = -Wall -O3 -o $(EXEC)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

main.o: $(SRCS)/main.cpp $(INCS)/weather_data.hpp $(INCS)/color_tables.hpp
	$(CC) $(CFLAGS) $(SRCS)/main.cpp

cell.o: $(SRCS)/cell.cpp $(INCS)/cell.hpp
	$(CC) $(CFLAGS) $(SRCS)/cell.cpp

weather_data.o: $(SRCS)/weather_data.cpp $(INCS)/weather_data.hpp $(LIBS)/lodepng.h
	$(CC) $(CFLAGS) $(SRCS)/weather_data.cpp

lodepng.o: $(LIBS)/lodepng.cpp $(LIBS)/lodepng.h
	$(CC) $(CFLAGS) $(LIBS)/lodepng.cpp

clean:
	@rm *.o