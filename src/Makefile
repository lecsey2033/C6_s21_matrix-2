CC=gcc
CFLAGS=-Wall -Wextra -Werror -g -std=c11
LIBS=-lcheck
LFLAGS=-fprofile-arcs -ftest-coverage

all: clean s21_matrix.a test

s21_matrix.a: s21_matrix.o
	ar rcs s21_matrix.a s21_matrix.o
	cp s21_matrix.a libs21_matrix.a
	cp s21_matrix.a ls21_matrix.a
	ranlib libs21_matrix.a
	ranlib ls21_matrix.a

s21_matrix.o: s21_matrix.c
	$(CC) $(CFLAGS) -c s21_matrix.c -o s21_matrix.o

test: s21_matrix.a
	$(CC) $(CFLAGS) test.c s21_matrix.a -o test.out $(LIBS)
	./test.out

gcov_report: s21_matrix.c
	$(CC) $(CFLAGS) $(LFLAGS) test.c s21_matrix.c -o test.out $(LIBS)
	./test.out
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info

open:
	open report/index.html

rebuild: clean all

check:
	cp ../materials/linters/CPPLINT.cfg ./
	cppcheck --enable=all --suppress=missingIncludeSystem *.c *.h
	python3 ../materials/linters/cpplint.py --extensions=c *.c *.h
	make test
	CK_FORK=no leaks --atExit -- ./test.out
	make clean
	rm -rf CPPLINT.cfg

clean:
	rm -rf *.gc* *.o *.a *.out *.info test.out
	rm -rf report