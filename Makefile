default: all

all:
	gcc -O3 -std=c11 -march=native -o sudoku sudoku.c -Wall

omp:
	gcc -O3 -std=c11 -march=native sudoku.c -o sudoku -fopenmp -Wall

mpi:
	mpicc -O3 -o sudoku.c sudoku
	mpirun -np 4 ./sudoku

omp_mpi:
	mpicc -O3 -std=c11 -march=native sudoku.c -o sudoku -fopenmp

omp_mpi_run:
	mpirun -np $(np) ./sudoku $(file)

clean:
	rm sudoku
