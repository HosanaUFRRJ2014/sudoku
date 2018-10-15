default: sequential

sequential:
	gcc -O3 -std=c11 -march=native -o sudoku sudoku.c -Wall

thread:
	gcc -O3 -std=c11 -march=native -o sudoku sudoku.c -lpthread -Wall


clean:
	rm sudoku
