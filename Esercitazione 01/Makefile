CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS =

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)

main01.exe: main01.o Funzioni.o random.o
	$(CC) $^ -o $@ $(LIBS)

main02.exe: main02.o Funzioni.o random.o
	$(CC) $^ -o $@ $(LIBS)

main03.exe: main03.o Funzioni.o random.o
	$(CC) $^ -o $@ $(LIBS)

clean :
	rm *.o main0*.exe seed.out
