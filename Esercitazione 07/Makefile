CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe

Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)

main7.1.exe: main7.1.o random.o
	$(CC) random.o main7.1.o -o main7.1.exe

main7.1.o : main7.1.cpp
	$(CC) -c main7.1.cpp -o main7.1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Monte_Carlo_NVT.exe main7.1.exe
