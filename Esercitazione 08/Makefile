CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS =

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)

main01.exe: main01.o random.o
	$(CC) $^ -o $@ $(LIBS)

clean :
	rm *.o main0*.exe seed.out
