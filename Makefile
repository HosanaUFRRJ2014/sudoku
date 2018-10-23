default: all

all:
	gcc -O3 -std=c11 -march=native -o sudoku sudoku.c -Wall

omp:
	gcc -O3 -std=c11 -march=native sudoku.c -o sudoku -fopenmp -Wall

clean:
	rm sudoku
