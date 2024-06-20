NAME := s21_matrix.a
CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11
TEST_LIBS = -lcheck -lm -lpthread
SRC := s21_matrix.c 
SRC_TEST := test.c
OBJ := $(patsubst %.c,%.o,$(SRC))
GCOVFLAGS = -fprofile-arcs -ftest-coverage
OS := $(shell uname)

ifeq ($(OS),Linux)
	CFLAGS+=-D LINUX
	TST_LIBS+= -lcheck_pic -lrt -lsubunit
else
	CFLAGS+=-D MAC
endif

.PHONY: all test gcov_report clean rebuild check

all: $(NAME)

$(NAME): $(OBJ)
	ar rcs $@ $^

%.o: %.c s21_matrix.h
	 $(CC) $(CFLAGS) -c $< -o $@

test: $(SRC_TEST) $(NAME)
	$(CC) $< $(NAME) $(TEST_LIBS) -o test
	./test

gcov_report: $(NAME)
	$(CC) $(CFLAGS) --coverage $(SRC) $(SRC_TEST) $(TEST_LIBS) -o gcov_test
	chmod +x *
	./gcov_test
	lcov -t "gcov_test" -o gcov_test.info --no-external -c -d .
	genhtml -o report/ gcov_test.info

clean:
	rm -rf ./*.o ./*.a ./a.out gcov_test *.html *.css ./GcovReport ./*.gcno ./*.gcda ./report ./*.info ./*.dSYM ./CPPLINT* test

rebuild: clean s21_matrix.a gcov_report
	open ./report/index.html

check: s21_matrix.a
	cppcheck *.h *.c
	mv ../materials/linters/.clang-format .clang-format
	clang-format -i *.c *.h
	clang-format -n *.c *.h
	mv .clang-format ../materials/linters/.clang-format
	$(CC) test.c $(NAME) -lcheck
	$(CC) $(CFLAGS) -c test.c -o test.o
	$(CC) test.o $(NAME) $(TEST_LIBS) -o test
	CK_FORK=no leaks --atExit -- ./test