CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS =

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)

main.exe: main.o random.o GeneticAlg.h
	$(CC) $^ -o $@ $(LIBS)

clean :
	rm *.o main.exe seed.out
