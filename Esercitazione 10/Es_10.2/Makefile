# scrivere nel terminale "module load mpi/mpich-3.2-x86_64" prima di lanciare il make
CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11
LIBS =

main.exe: main.o random.o
	$(CC) $^ -o $@ $(LIBS)

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS)
	
main.o : main.cpp GeneticAlg.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

exec: main.exe
	mpiexec -np 4 ./main.exe file.txt
# se non metto "file.txt" il comando mi dà errore

clean :
	rm *.o main.exe seed.out
