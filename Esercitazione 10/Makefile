CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS =

main.exe: main.o random.o
	$(CC) $^ -o $@ $(LIBS)

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)

main.o : main.cpp GeneticAlg.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
